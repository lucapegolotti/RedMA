#include "BuildingBlock.hpp"
#include <iostream>

namespace RedMA
{

GeometricFace::
GeometricFace() :
  M_center(Vector3D(0.0,0.0,0.0)),
  M_normal(Vector3D(1.0,0.0,0.0)),
  M_radius(1.0),
  M_flag(0),
  M_diskFlag(0)
{
}

GeometricFace::
GeometricFace(Vector3D center, Vector3D normal, double radius,
              unsigned int flag) :
  M_center(center),
  M_normal(normal),
  M_radius(radius),
  M_flag(flag),
  M_diskFlag(0)
{
}

GeometricFace::
GeometricFace(Vector3D center, Vector3D normal, double radius,
              unsigned int flag, unsigned int flagDisk) :
  M_center(center),
  M_normal(normal),
  M_radius(radius),
  M_flag(flag),
  M_diskFlag(flagDisk)
{
}

void
BuildingBlock::
setDatafile(const GetPot& datafile)
{
    M_datafile = datafile;
}

void
GeometricFace::
print() const
{
    printlog(WHITE, "[GeometricFace]\n");
    printlog(WHITE, std::string("\tcenter = (") + std::to_string(M_center[0]) +
                    "," + std::to_string(M_center[1]) + "," +
                    std::to_string(M_center[2]) + ")\n");
    printlog(WHITE, std::string("\tnormal = (") + std::to_string(M_normal[0]) +
                    "," + std::to_string(M_normal[1]) + "," +
                    std::to_string(M_normal[2]) + ")\n");
    printlog(WHITE, std::string("\tradius = ") + std::to_string(M_radius) + "\n");
    printlog(WHITE, std::string("\tflag = ") + std::to_string(M_flag) + "\n");

}

BuildingBlock::
BuildingBlock(commPtr_Type comm, std::string refinement, bool verbose) :
  M_comm(comm),
  M_refinement(refinement),
  M_verbose(verbose),
  M_isChild(false)
{
    if (M_comm->MyPID() != 0)
        M_verbose = false;

    double infty = GeometricParametersHandler::infty;
    double mp2 = 20 * M_PI;

    // rotation axis and angle
    M_parametersHandler.registerParameter("rotation_axis_x", 1.0, -infty, infty);
    M_parametersHandler.registerParameter("rotation_axis_y", 0.0, -infty, infty);
    M_parametersHandler.registerParameter("rotation_axis_z", 0.0, -infty, infty);
    M_parametersHandler.registerParameter("alpha", 0.0, -mp2, mp2, false, true);


    M_parametersHandler.registerParameter("alpha_axis", 0.0, -mp2, mp2, true, true);

    // scale
    M_parametersHandler.registerParameter("scale", 1.0, 0.0, infty);

    // translation
    M_parametersHandler.registerParameter("bx", 0.0, -infty, infty);
    M_parametersHandler.registerParameter("by", 0.0, -infty, infty);
    M_parametersHandler.registerParameter("bz", 0.0, -infty, infty);
}

int
BuildingBlock::
setParameterValue(std::string key, double value)
{
    return M_parametersHandler.setParameterValue(key, value);
}

std::map<std::string,std::shared_ptr<GeometricParameter> >&
BuildingBlock::
getParametersMap()
{
    return M_parametersHandler.getParametersMap();
}

int
BuildingBlock::
readMesh(std::string meshdir)
{
    printlog(MAGENTA, "[" + M_name +
                    " BuildingBlock] reading mesh ...\n",
                    M_verbose);

    meshPtr_Type fullMesh(new mesh_Type(M_comm));
    LifeV::MeshData meshData;
    GetPot meshDatafile(meshdir + M_datafileName);
    meshDatafile.set("mesh/mesh_file", M_meshName.c_str());
    meshData.setup(meshDatafile, "mesh");
    meshData.setMeshDir(meshdir);
    LifeV::readMesh(*fullMesh,meshData);

    LifeV::MeshPartitioner<mesh_Type> meshPart(fullMesh, M_comm);

    // small trick to redirect std cout
    // CoutRedirecter ct;
    // ct.redirect();
    // meshPart.doPartition(fullMesh, M_comm);
    M_mesh.reset(new mesh_Type(M_comm));
    M_mesh = meshPart.meshPartition();

    // printlog(CYAN, ct.restore(), M_verbose);
    printlog(MAGENTA, "done\n", M_verbose);

    return 0;
}

std::string
BuildingBlock::
name()
{
    return M_name;
}

BuildingBlock::Matrix3D
BuildingBlock::
computeRotationMatrix(unsigned int axis, double angle)
{
    Matrix3D R;
    if (axis == 0)
    {
        R(0,0) = 1.;
        R(0,1) = 0.;
        R(0,2) = 0.;
        R(1,0) = 0.;
        R(1,1) = std::cos(angle);
        R(1,2) = -std::sin(angle);
        R(2,0) = 0.;
        R(2,1) = std::sin(angle);
        R(2,2) = std::cos(angle);
    }
    else if (axis == 1)
    {
        R(0,0) = std::cos(angle);
        R(0,1) = 0.;
        R(0,2) = std::sin(angle);
        R(1,0) = 0.;
        R(1,1) = 1.;
        R(1,2) = 0.;
        R(2,0) = -std::sin(angle);
        R(2,1) = 0.;
        R(2,2) = std::cos(angle);
    }
    else if (axis == 2)
    {
        R(0,0) = std::cos(angle);
        R(0,1) = -std::sin(angle);
        R(0,2) = 0.;
        R(1,0) = std::sin(angle);
        R(1,1) = std::cos(angle);
        R(1,2) = 0.;
        R(2,0) = 0.;
        R(2,1) = 0.;
        R(2,2) = 1.;
    }
    return R;
}

BuildingBlock::Matrix3D
BuildingBlock::
computeRotationMatrix(Vector3D axis, double angle)
{
    if (axis.norm() < 1e-15 && std::abs(angle) < 1e-15)
    {
        axis[0] = 1;
        axis[1] = 0;
        axis[2] = 0;
    }
    else if (axis.norm() < 0)
    {
        throw new Exception("Rotation axis is set to (0,0,0)!");
    }

    axis = axis / axis.norm();

    Matrix3D R;
    double mcos = std::cos(angle);
    double omcos = 1.0 - mcos;
    double msin = std::sin(angle);

    // enforce normality of axis
    axis = axis / axis.norm();

    R(0,0) = mcos + axis[0] * axis[0] * omcos;
    R(0,1) = axis[0] * axis[1] * omcos - axis[2] * msin;
    R(0,2) = axis[0] * axis[2] * omcos + axis[1] * msin;
    R(1,0) = axis[0] * axis[1] * omcos + axis[2] * msin;
    R(1,1) = mcos + axis[1] * axis[1] * omcos;
    R(1,2) = axis[1] * axis[2] * omcos - axis[0] * msin;
    R(2,0) = axis[2] * axis[0] * omcos - axis[1] * msin;
    R(2,1) = axis[1] * axis[2] * omcos + axis[0] * msin;
    R(2,2) = mcos + axis[2] * axis[2] * omcos;

    return R;
}

void
BuildingBlock::
applyAffineTransformation(bool transformMesh)
{
    printlog(MAGENTA, "[" + M_name +
                    " BuildingBlock] applying affine transformation ...\n",
                    M_verbose);
    if (transformMesh && !M_mesh)
    {
        std::string errorMsg = std::string("[") + M_name + " BuildingBlock] " +
                               " mesh has not being read!";
        throw Exception(errorMsg);
    }

    std::shared_ptr<LifeV::MeshUtility::MeshTransformer<mesh_Type> > transformer;

    if (transformMesh)
        transformer.reset(new Transformer(*M_mesh));

    double scale;
    Vector3D scaleVec;
    Vector3D rotation;
    Vector3D translation;

    if (M_isChild)
    {
        printlog(GREEN, "[" + M_name +
                     " BuildingBlock] is child: conforming inlet/parent_outlet ...\n",
                     M_verbose);

        M_R =  computeRotationMatrix(M_inletRotationAxis, M_inletAngle);

        translation = M_inletTranslation;
        scale = M_inletScale;

        auto foo = std::bind(rotationFunction,
                             std::placeholders::_1,
                             std::placeholders::_2,
                             std::placeholders::_3,
                             M_R, translation, scale);

        std::string axisStr = std::string("(");
        axisStr = axisStr + std::to_string(M_inletRotationAxis[0]) + ",";
        axisStr = axisStr + std::to_string(M_inletRotationAxis[1]) + ",";
        axisStr = axisStr + std::to_string(M_inletRotationAxis[2]) + ")";

        printlog(YELLOW, "[" + M_name +
                     " BuildingBlock] rotating about axis " + axisStr + " and" +
                     " angle " + std::to_string(M_inletAngle) +
                     " ...\n", M_verbose);

        std::string transStr = std::string("(");
        transStr  = "(";
        transStr = transStr + std::to_string(translation[0]) + ",";
        transStr = transStr + std::to_string(translation[1]) + ",";
        transStr = transStr + std::to_string(translation[2]) + ")";

        printlog(YELLOW, "[" + M_name +
                  " BuildingBlock] translating with vector " + transStr +
                  " ...\n", M_verbose);

        printlog(YELLOW, "[" + M_name +
                " BuildingBlock] applying scaling of " + std::to_string(scale) +
                " ...\n", M_verbose);

        if (transformMesh)
            transformer->transformMesh(foo);
        printlog(GREEN, "done\n", M_verbose);
    }
    else
    {
        printlog(GREEN, "[" + M_name +
                     " BuildingBlock] is root: applying initial affine transformation ...\n",
                     M_verbose);
        scale = M_parametersHandler["scale"];

        rotation[0] = M_parametersHandler["rotation_axis_x"];
        rotation[1] = M_parametersHandler["rotation_axis_y"];
        rotation[2] = M_parametersHandler["rotation_axis_z"];

        double angle = M_parametersHandler["alpha"];

        translation[0] = M_parametersHandler["bx"];
        translation[1] = M_parametersHandler["by"];
        translation[2] = M_parametersHandler["bz"];

        std::string axisStr = std::string("(");
        axisStr = axisStr + std::to_string(rotation[0]) + ",";
        axisStr = axisStr + std::to_string(rotation[1]) + ",";
        axisStr = axisStr + std::to_string(rotation[2]) + ")";

        printlog(YELLOW, "[" + M_name +
                     " BuildingBlock] rotating about axis " + axisStr +
                     " with angle " + std::to_string(angle) + " ...\n", M_verbose);


        std::string transStr = std::string("(");
        transStr  = "(";
        transStr = transStr + std::to_string(translation[0]) + ",";
        transStr = transStr + std::to_string(translation[1]) + ",";
        transStr = transStr + std::to_string(translation[2]) + ")";

        printlog(YELLOW, "[" + M_name +
               " BuildingBlock] translating with vector " + transStr +
               " ...\n", M_verbose);

        M_R =  computeRotationMatrix(rotation, angle);

        auto foo = std::bind(rotationFunction,
                             std::placeholders::_1,
                             std::placeholders::_2,
                             std::placeholders::_3,
                             M_R, translation, scale);

        if (transformMesh)
            transformer->transformMesh(foo);

        printlog(YELLOW, "[" + M_name +
                " BuildingBlock] applying scaling of " + std::to_string(scale) +
                " ...\n", M_verbose);

        printlog(GREEN, "done\n", M_verbose);
    }

    applyAffineTransformationGeometricFace(M_inlet,M_R,translation,scale);
    for (std::vector<GeometricFace>::iterator it = M_outlets.begin();
         it != M_outlets.end(); it++)
    {
        applyAffineTransformationGeometricFace(*it, M_R, translation,scale);
    }

    // Handle rotation along the axis of the inlet
    double angle = M_parametersHandler["alpha_axis"];

    Matrix3D Raxis = computeRotationMatrix(M_inlet.M_normal, angle);
    Vector3D transZero = M_inlet.M_center - Raxis * M_inlet.M_center;

    auto foo = std::bind(rotationFunction,
                         std::placeholders::_1,
                         std::placeholders::_2,
                         std::placeholders::_3,
                         Raxis, transZero, 1.0);

    if (transformMesh)
        transformer->transformMesh(foo);

    for (std::vector<GeometricFace>::iterator it = M_outlets.begin();
         it != M_outlets.end(); it++)
    {
        applyAffineTransformationGeometricFace(*it, Raxis, transZero, 1.0);
    }

    printlog(MAGENTA, "done\n", M_verbose);
}

void
BuildingBlock::
applyAffineTransformationGeometricFace(GeometricFace& face,
                                       const Matrix3D& affineMatrix,
                                       const Vector3D& translation,
                                       const double& scale)
{
    face.M_center = scale * affineMatrix * face.M_center;
    face.M_center = face.M_center + translation;

    face.M_normal = affineMatrix * face.M_normal;

    face.M_radius = face.M_radius * scale;
}

void
BuildingBlock::
rotationFunction(double& x, double& y, double& z, const Matrix3D& affMatrix,
               const Vector3D& transl, const double& scale)
{
    Vector3D vec;
    vec[0] = x; vec[1] = y; vec[2] = z;

    vec = vec * scale;
    vec = affMatrix * vec + transl;

    x = vec[0];
    y = vec[1];
    z = vec[2];
}

void
BuildingBlock::
computeRotationAxisAndAngle(Vector3D vectorToMove,
                            Vector3D vectorToReach,
                            Vector3D& axis,
                            double& alpha)
{

    axis = vectorToMove.cross(vectorToReach);
    axis = axis / axis.norm();

    alpha = std::acos(vectorToMove.dot(vectorToReach) /
                      (vectorToMove.norm() * vectorToReach.norm()));
}

void
BuildingBlock::
dumpMesh(std::string outdir, std::string meshdir, std::string outputName)
{
    printlog(MAGENTA, "[" + M_name +
                    " BuildingBlock] dumping mesh to file ...\n",
                    M_verbose);

    if (!M_mesh)
    {
        std::string msg = "Mesh has not been read yet!\n";
        throw Exception(msg);
    }
    boost::filesystem::create_directory(outdir);

    GetPot exporterDatafile(meshdir + M_datafileName);
    LifeV::ExporterVTK<mesh_Type> exporter(exporterDatafile, outputName);
    exporter.setMeshProcId(M_mesh, M_comm->MyPID());

    FESpacePtr_Type dummyFespace(new FESpace_Type(M_mesh, "P1", 3, M_comm));
    vectorPtr_Type zero(new vector_Type(dummyFespace->map()) );
    zero->zero();

    exporter.addVariable(LifeV::ExporterData<mesh_Type>::ScalarField, "z",
                         dummyFespace, zero, 0);
    exporter.setPostDir(outdir);

    CoutRedirecter ct;
    ct.redirect();
    exporter.postProcess(0.0);
    printlog(CYAN, ct.restore(), M_verbose);

    printlog(MAGENTA, "done\n", M_verbose);
}

GeometricFace
BuildingBlock::
getOutlet(unsigned int outletIndex) const
{
    if (outletIndex >= M_outlets.size())
    {
        std::string msg = "Requesting access to outlet that does not exist!";
        throw Exception(msg);
    }

    return M_outlets[outletIndex];
}

std::vector<GeometricFace>
BuildingBlock::
getOutlets() const
{
    return M_outlets;
}

GeometricFace
BuildingBlock::
getInlet() const
{
    return M_inlet;
}

void
BuildingBlock::
mapChildInletToParentOutlet(GeometricFace parentOutlet)
{
    M_isChild = true;

    Vector3D iNormal = -1 * M_inlet.M_normal;
    Vector3D oNormal = parentOutlet.M_normal;
    Vector3D iCenter = M_inlet.M_center;
    Vector3D oCenter = parentOutlet.M_center;

    M_inletScale = parentOutlet.M_radius / M_inlet.M_radius;

    M_inletTranslation = oCenter - iCenter;

    M_inletRotationAxis = iNormal.cross(oNormal);
    // in this case we don't have to rotate
    if (M_inletRotationAxis.norm() < 1e-15)
    {
        M_inletRotationAxis = oNormal;
        M_inletAngle = 0;
    }
    else
    {
        M_inletRotationAxis = M_inletRotationAxis / M_inletRotationAxis.norm();
        M_inletAngle = std::acos(iNormal.dot(oNormal) /
                                (iNormal.norm() * oNormal.norm()));
    }
}

void
BuildingBlock::
applyGlobalTransformation(bool transformMesh)
{
    applyNonAffineTransformation(transformMesh);
    applyAffineTransformation(transformMesh);
}

void
BuildingBlock::
setIsChild(bool isChild)
{
    M_isChild = isChild;
}

BuildingBlock::meshPtr_Type
BuildingBlock::
getMesh()
{
    return M_mesh;
}

void
BuildingBlock::
setRandom()
{
    M_parametersHandler.randomizeParameters();
}

double
BuildingBlock::
fZero(const double& t, const double& x, const double& y,
      const double& z, const LifeV::ID& i)
{
    return 0.0;
}

BuildingBlock::Matrix3D
BuildingBlock::
computeJacobianGlobalTransformation(const double& x,
                                    const double& y,
                                    const double& z)
{
    Matrix3D nonLinearJacobian = computeJacobianNonAffineTransformation(x,y,z);

    return M_inletScale * M_R * nonLinearJacobian;
}

}  // namespace BuildingBlock
