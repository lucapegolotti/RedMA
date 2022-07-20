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
  M_ringFlag(0)
{
}

GeometricFace::
GeometricFace(Vector3D center, Vector3D normal, double radius,
              unsigned int flag) :
  M_center(center),
  M_normal(normal),
  M_radius(radius),
  M_flag(flag),
  M_ringFlag(0)
{
}

GeometricFace::
GeometricFace(Vector3D center, Vector3D normal, double radius,
              unsigned int flag, unsigned int flagRing) :
  M_center(center),
  M_normal(normal),
  M_radius(radius),
  M_flag(flag),
  M_ringFlag(flagRing)
{
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
    printlog(WHITE, std::string("\tring flag = ") + std::to_string(M_ringFlag) + "\n");

}

bool
operator==(const GeometricFace& lhs, const GeometricFace& rhs)
{
    typedef LifeV::VectorSmall<3>   Vector3D;

    Vector3D normalLhs;
    Vector3D centerLhs;
    Vector3D normalRhs;
    Vector3D centerRhs;

    normalLhs = lhs.M_normal;
    normalLhs *= (-1);
    centerLhs = lhs.M_center;

    normalRhs = rhs.M_normal;
    centerRhs = rhs.M_center;

    double diff = (normalLhs - normalRhs).norm() + (centerLhs - centerRhs).norm();

    return (diff <= 1e-3);
}

BuildingBlock::
BuildingBlock(EPETRACOMM comm, std::string refinement, bool verbose) :
  M_comm(comm),
  M_refinement(refinement),
  M_verbose(verbose),
  M_isChild(false),
  M_discrMethod("none"),
  M_assemblerType("none"),
  M_wallFlag(10)
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

    // angle around axis
    M_parametersHandler.registerParameter("alpha_axis", 0.0, -mp2, mp2, true, true);

    // scale
    M_parametersHandler.registerParameter("scale", 1.0, 0.0, infty);

    // translation
    M_parametersHandler.registerParameter("bx", 0.0, -infty, infty);
    M_parametersHandler.registerParameter("by", 0.0, -infty, infty);
    M_parametersHandler.registerParameter("bz", 0.0, -infty, infty);

    // M_displacement->zero();
}

void
BuildingBlock::
setDatafile(const DataContainer& datafile)
{
    M_datafile = datafile;
}

int
BuildingBlock::
setParameterValue(std::string key, double value)
{
    return M_parametersHandler.setParameterValue(key, value);
}

void
BuildingBlock::
matrixInverse(Matrix3D& matrix, Matrix3D* inverse)
{
    double& a = matrix(0,0);
    double& b = matrix(0,1);
    double& c = matrix(0,2);
    double& d = matrix(1,0);
    double& e = matrix(1,1);
    double& f = matrix(1,2);
    double& g = matrix(2,0);
    double& h = matrix(2,1);
    double& i = matrix(2,2);

    (*inverse)(0,0) = (e*i - f*h)/(a*e*i - a*f*h - b*d*i + b*f*g + c*d*h - c*e*g);
    (*inverse)(0,1) = -(b*i - c*h)/(a*e*i - a*f*h - b*d*i + b*f*g + c*d*h - c*e*g);
    (*inverse)(0,2) = (b*f - c*e)/(a*e*i - a*f*h - b*d*i + b*f*g + c*d*h - c*e*g);
    (*inverse)(1,0) = -(d*i - f*g)/(a*e*i - a*f*h - b*d*i + b*f*g + c*d*h - c*e*g);
    (*inverse)(1,1) = (a*i - c*g)/(a*e*i - a*f*h - b*d*i + b*f*g + c*d*h - c*e*g);
    (*inverse)(1,2) = -(a*f - c*d)/(a*e*i - a*f*h - b*d*i + b*f*g + c*d*h - c*e*g);
    (*inverse)(2,0) = (d*h - e*g)/(a*e*i - a*f*h - b*d*i + b*f*g + c*d*h - c*e*g);
    (*inverse)(2,1) = -(a*h - b*g)/(a*e*i - a*f*h - b*d*i + b*f*g + c*d*h - c*e*g);
    (*inverse)(2,2) = (a*e - b*d)/(a*e*i - a*f*h - b*d*i + b*f*g + c*d*h - c*e*g);
}

std::map<std::string,shp<GeometricParameter> >&
BuildingBlock::
getParametersMap()
{
    return M_parametersHandler.getParametersMap();
}

void
BuildingBlock::
readMesh(std::string meshdir)
{
    printlog(MAGENTA, "[" + M_name +
                    " BuildingBlock] reading mesh ...\n",
                    M_verbose);

    shp<MESH> fullMesh(new MESH(M_comm));
    LifeV::MeshData meshData;
    GetPot meshDatafile(meshdir + M_datafileName);
    meshDatafile.set("mesh/mesh_file", M_meshName.c_str());
    meshData.setup(meshDatafile, "mesh");
    meshData.setMeshDir(meshdir);
    LifeV::readMesh(*fullMesh,meshData);

    LifeV::MeshPartitioner<MESH> meshPart(fullMesh, M_comm);

    M_mesh.reset(new MESH(M_comm));
    M_mesh = meshPart.meshPartition();

    printlog(MAGENTA, "done\n", M_verbose);
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

BuildingBlock::Matrix3D
BuildingBlock::
computeLocalRotationMatrix(Vector3D vec, const unsigned int& index) const
{
    if (std::abs(vec.norm() - 1.0) <= 1e-10)
        vec.normalize();

    if (index > 2)
        printlog(YELLOW, "[computeRotationMatrix] WARNING: "
                         "the second argument must be either 0, 1 or 2!"
                         "Setting it to 2 as default value",
                 this->M_datafile.getVerbose());

    Matrix3D mat;

    std::pair<Vector3D, Vector3D> tangents = this->computeLocalTangentVersors(vec);

    // assign the prescribed normalized vector to the prescribed column
    mat[0][index] = vec[0];
    mat[1][index] = vec[1];
    mat[2][index] = vec[2];

    if (index != 0) {
        mat[0][0] = tangents.first[0];
        mat[1][0] = tangents.first[1];
        mat[2][0] = tangents.first[2];
    }
    else {
        mat[0][1] = tangents.first[0];
        mat[1][1] = tangents.first[1];
        mat[2][1] = tangents.first[2];
    }

    if (index != 2) {
        mat[0][2] = tangents.second[0];
        mat[1][2] = tangents.second[1];
        mat[2][2] = tangents.second[2];
    }
    else {
        mat[0][1] = tangents.second[0];
        mat[1][1] = tangents.second[1];
        mat[2][1] = tangents.second[2];
    }

    return mat;
}

std::pair<BuildingBlock::Vector3D, BuildingBlock::Vector3D>
BuildingBlock::
computeLocalTangentVersors(const Vector3D& normal)
{
    Vector3D t1;
    Vector3D t2;
    std::pair<Vector3D, Vector3D> res;

    double nx = normal[0];
    double ny = normal[1];
    double nz = normal[2];

    double nx2 = std::sqrt(ny * ny + nz * nz);
    double ny2 = std::sqrt(nx * nx + nz * nz);
    double nz2 = std::sqrt(nx * nx + ny * ny);

    if ((nx2 >= ny2) && (nx2 >= nz2))
    {
        //We create t1
        t1[0] = 0;
        t1[1] = nz / nx2;
        t1[2] = -ny / nx2;

        //We create t2
        t2[0] = -nx2;
        t2[1] = nx * ny / nx2;
        t2[2] = nx * nz / nx2;
    }
    else if ((ny2 >= nx2) && (ny2 >= nz2))
    {
        //We create t1
        t1[0] = -nz / ny2;
        t1[1] = 0;
        t1[2] = nx / ny2;

        //We create t2
        t2[0] = nx * ny / ny2;
        t2[1] = -ny2;
        t2[2] = ny * nz / ny2;
    }
    else
    {
        //We create t1
        t1[0] = ny / nz2;
        t1[1] = -nx / nz2;
        t1[2] = 0;

        //We create t2
        t2[0] = nx * nz / nz2;
        t2[1] = ny * nz / nz2;
        t2[2] = -nz2;
    }

    res.first = t1;
    res.second = t2;

    return res;
}

void
BuildingBlock::
computeMembraneThickness() {

    M_membraneThicknessComputer.reset(new MembraneThicknessComputer(M_datafile,
                                                                    M_comm));

    std::vector<double> radia_in;
    std::vector<unsigned int> flags_in;
    for (auto in : M_inlets)
    {
        radia_in.push_back(in.M_radius);
        flags_in.push_back(in.M_ringFlag);
    }

    std::vector<double> radia_out;
    std::vector<unsigned int> flags_out;
    for (auto out : M_outlets)
    {
        radia_out.push_back(out.M_radius);
        flags_out.push_back(out.M_ringFlag);
    }

    M_membraneThicknessComputer->setup(radia_in, radia_out,
                                       flags_in, flags_out,
                                       M_wallFlag, M_mesh);

    M_membraneThicknessComputer->solve();
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

    shp<LifeV::MeshUtility::MeshTransformer<MESH> > transformer;

    if (transformMesh)
        transformer.reset(new Transformer(*M_mesh));

    Vector3D rotation;

    if (M_isChild)
    {
        printlog(GREEN, "[" + M_name +
                     " BuildingBlock] is child: conforming inlet/parent_outlet ...\n",
                     M_verbose);

        M_R =  computeRotationMatrix(M_inletRotationAxis, M_inletAngle);

        M_translation = M_inletTranslation;
        M_scale = M_inletScale;

        auto foo = std::bind(rotationFunction,
                             std::placeholders::_1,
                             std::placeholders::_2,
                             std::placeholders::_3,
                             M_R, M_translation, M_scale);

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
        transStr = transStr + std::to_string(M_translation[0]) + ",";
        transStr = transStr + std::to_string(M_translation[1]) + ",";
        transStr = transStr + std::to_string(M_translation[2]) + ")";

        printlog(YELLOW, "[" + M_name +
                  " BuildingBlock] translating with vector " + transStr +
                  " ...\n", M_verbose);

        printlog(YELLOW, "[" + M_name +
                " BuildingBlock] applying scaling of " + std::to_string(M_scale) +
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
        M_scale = M_parametersHandler["scale"];

        rotation[0] = M_parametersHandler["rotation_axis_x"];
        rotation[1] = M_parametersHandler["rotation_axis_y"];
        rotation[2] = M_parametersHandler["rotation_axis_z"];

        double angle = M_parametersHandler["alpha"];

        M_translation[0] = M_parametersHandler["bx"];
        M_translation[1] = M_parametersHandler["by"];
        M_translation[2] = M_parametersHandler["bz"];

        std::string axisStr = std::string("(");
        axisStr = axisStr + std::to_string(rotation[0]) + ",";
        axisStr = axisStr + std::to_string(rotation[1]) + ",";
        axisStr = axisStr + std::to_string(rotation[2]) + ")";

        printlog(YELLOW, "[" + M_name +
                     " BuildingBlock] rotating about axis " + axisStr +
                     " with angle " + std::to_string(angle) + " ...\n", M_verbose);


        std::string transStr = std::string("(");
        transStr  = "(";
        transStr = transStr + std::to_string(M_translation[0]) + ",";
        transStr = transStr + std::to_string(M_translation[1]) + ",";
        transStr = transStr + std::to_string(M_translation[2]) + ")";

        printlog(YELLOW, "[" + M_name +
               " BuildingBlock] translating with vector " + transStr +
               " ...\n", M_verbose);

        printlog(YELLOW, "[" + M_name +
                         " BuildingBlock] applying scaling of " + std::to_string(M_scale) +
                         " ...\n", M_verbose);

        M_R =  computeRotationMatrix(rotation, angle);

        auto foo = std::bind(rotationFunction,
                             std::placeholders::_1,
                             std::placeholders::_2,
                             std::placeholders::_3,
                             M_R, M_translation, M_scale);

        if (transformMesh)
            transformer->transformMesh(foo);

        printlog(GREEN, "done\n", M_verbose);
    }

    for (std::vector<GeometricFace>::iterator it = M_inlets.begin();
         it != M_inlets.end(); it++)
    {
        applyAffineTransformationGeometricFace(*it, M_R, M_translation, M_scale);
    }
    for (std::vector<GeometricFace>::iterator it = M_outlets.begin();
         it != M_outlets.end(); it++)
    {
        applyAffineTransformationGeometricFace(*it, M_R, M_translation, M_scale);
    }

    // Handle rotation along the axis of the inlet 0
    double angle = M_parametersHandler["alpha_axis"];

    M_Raxis = computeRotationMatrix(M_inlets[0].M_normal, angle);
    Vector3D transZero = M_inlets[0].M_center - M_Raxis * M_inlets[0].M_center;

    auto foo = std::bind(rotationFunction,
                         std::placeholders::_1,
                         std::placeholders::_2,
                         std::placeholders::_3,
                         M_Raxis, transZero, 1.0);

    if (transformMesh)
        transformer->transformMesh(foo);

    for (std::vector<GeometricFace>::iterator it = M_inlets.begin();
         it != M_inlets.end(); it++)
    {
        applyAffineTransformationGeometricFace(*it, M_Raxis, transZero, 1.0);
    }

    for (std::vector<GeometricFace>::iterator it = M_outlets.begin();
         it != M_outlets.end(); it++)
    {
        applyAffineTransformationGeometricFace(*it, M_Raxis, transZero, 1.0);
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
    fs::create_directory(outdir);

    GetPot exporterDatafile(meshdir + M_datafileName);
    LifeV::ExporterVTK<MESH> exporter(exporterDatafile, outputName);
    exporter.setMeshProcId(M_mesh, M_comm->MyPID());

    shp<FESPACE> dummyFespace(new FESPACE(M_mesh, "P1", 3, M_comm));
    shp<VECTOREPETRA> zero(new VECTOREPETRA(dummyFespace->map()) );
    zero->zero();

    exporter.addVariable(LifeV::ExporterData<MESH>::ScalarField, "z",
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
        throw new Exception(msg);
    }

    return M_outlets[outletIndex];
}

GeometricFace
BuildingBlock::
getInlet(unsigned int inletIndex) const
{
    if (inletIndex >= M_inlets.size())
    {
        std::string msg = "Requesting access to inlet that does not exist!";
        throw new Exception(msg);
    }

    return M_inlets[inletIndex];
}

std::vector<GeometricFace>
BuildingBlock::
getOutlets() const
{
    return M_outlets;
}

std::vector<GeometricFace>
BuildingBlock::
getInlets() const
{
    return M_inlets;
}

void
BuildingBlock::
mapChildInletToParentOutlet(GeometricFace parentOutlet)
{
    M_isChild = true;

    Vector3D iNormal = -1 * M_inlets[0].M_normal;
    Vector3D oNormal = parentOutlet.M_normal;
    Vector3D iCenter = M_inlets[0].M_center;
    Vector3D oCenter = parentOutlet.M_center;

    M_inletScale = parentOutlet.M_radius / M_inlets[0].M_radius;

    M_inletTranslation = oCenter - iCenter;

    M_inletRotationAxis = iNormal.cross(oNormal);
    // in this case we don't have to rotate
    if (M_inletRotationAxis.norm() < 1e-5)
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

shp<MESH>
BuildingBlock::
getMesh() const
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
    return M_scale * M_Raxis * M_R * nonLinearJacobian;
}

void
BuildingBlock::
globalTransf(double& x, double& y, double& z)
{
    nonAffineTransf(x,y,z);
    affineTransf(x,y,z);
}

void
BuildingBlock::
affineTransf(double& x, double& y, double& z)
{
    Vector3D coor(x,y,z);
    Vector3D res = M_scale * M_Raxis * M_R * coor + M_translation;

    x = res(0); y = res(1); z = res(2);
}

}  // namespace BuildingBlock
