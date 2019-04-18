#include "BuildingBlock.hpp"
#include <iostream>

namespace RedMA
{

GeometricFace::
GeometricFace() :
  M_center(Vector3D(0.0,0.0,0.0)),
  M_normal(Vector3D(1.0,0.0,0.0)),
  M_radius(1.0)
{
}

GeometricFace::
GeometricFace(Vector3D center, Vector3D normal, double radius) :
  M_center(center),
  M_normal(normal),
  M_radius(radius)
{
}

void
GeometricFace::
print()
{
    printlog(WHITE, "[GeometricFace]\n");
    printlog(WHITE, std::string("\tcenter = (") + std::to_string(M_center[0]) +
                    "," + std::to_string(M_center[1]) + "," +
                    std::to_string(M_center[2]) + ")\n");
    printlog(WHITE, std::string("\tnormal = (") + std::to_string(M_normal[0]) +
                    "," + std::to_string(M_normal[1]) + "," +
                    std::to_string(M_normal[2]) + ")\n");
    printlog(WHITE, std::string("\tradius = ") + std::to_string(M_radius) + "\n");
}

BuildingBlock::
BuildingBlock(commPtr_Type comm, bool verbose) :
  M_comm(comm),
  M_verbose(verbose),
  M_isChild(false)
{
    if (M_comm->MyPID() != 0)
        M_verbose = false;

    // rotation angle
    M_parametersMap["alphax"] = 0.0;
    M_parametersMap["alphay"] = 0.0;
    M_parametersMap["alphaz"] = 0.0;
    M_parametersMap["alpha_axis"] = 0.0;

    // scale
    M_parametersMap["scale"] = 1.0;

    // translation
    M_parametersMap["bx"] = 0.0;
    M_parametersMap["by"] = 0.0;
    M_parametersMap["bz"] = 0.0;
}

void
BuildingBlock::
setParameterValue(std::string key, double value)
{
    if (M_parametersMap.find(key) != M_parametersMap.end())
    {
        M_parametersMap[key] = value;
    }
    else
    {
        std::string errorMsg =
                     "Parameter with key " + key + " and value " +
                     std::to_string(value) + " is not contained" +
                     " in " + M_name + " building block!\n";

        throw Exception(errorMsg);
    }
}

std::map<std::string,double>&
BuildingBlock::
getParametersMap()
{
    return M_parametersMap;
}

int
BuildingBlock::
readMesh(std::string meshdir)
{
    printlog(GREEN, "[" + M_name +
                    " BuildingBlock] reading mesh ...\n",
                    M_verbose);

    meshPtr_Type fullMesh(new mesh_Type(M_comm));
    LifeV::MeshData meshData;
    GetPot meshDatafile(meshdir + "datafiles/" + M_datafileName);
    meshData.setup(meshDatafile, "mesh");
    meshData.setMeshDir(meshdir);

    LifeV::readMesh(*fullMesh,meshData);

    LifeV::MeshPartitioner<mesh_Type> meshPart;

    // small trick to redirect std cout
    CoutRedirecter ct;
    ct.redirect();
    meshPart.doPartition(fullMesh, M_comm);
    M_mesh = meshPart.meshPartition();

    printlog(CYAN, ct.restore(), M_verbose);
    printlog(GREEN, "done\n", M_verbose);

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
    Matrix3D R;
    double mcos = std::cos(angle);
    double omcos = 1.0 - mcos;
    double msin = std::sin(angle);

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
applyAffineTransformation()
{
    LifeV::MeshUtility::MeshTransformer<mesh_Type> transformer(*M_mesh);

    Matrix3D R, R1, R2, R3;
    double scale;
    Vector3D scaleVec;
    Vector3D rotation;
    Vector3D translation;

    if (M_isChild)
    {
        R =  computeRotationMatrix(M_inletRotationAxis, M_inletAngle);

        translation = M_inletTranslation;
        scale = M_inletScale;

        std::cout << "rotation axis" << std::endl;
        std::cout << M_inletRotationAxis[0] << std::endl;
        std::cout << M_inletRotationAxis[1] << std::endl;
        std::cout << M_inletRotationAxis[2] << std::endl;

        std::cout << "rotation angle" << std::endl;
        std::cout << M_inletAngle << std::endl;


        auto foo = std::bind(rotationFunction,
                             std::placeholders::_1,
                             std::placeholders::_2,
                             std::placeholders::_3,
                             R, translation, scale);

        transformer.transformMesh(foo);
    }
    else
    {
        scale = M_parametersMap["scale"];
        scaleVec[0] = scale; scaleVec[1] = scale; scaleVec[2] = scale;

        // with minus in front: counterclockwise rotation
        rotation[0] = -M_parametersMap["alphax"];
        rotation[1] = -M_parametersMap["alphay"];
        rotation[2] = -M_parametersMap["alphaz"];

        translation[0] = M_parametersMap["bx"];
        translation[1] = M_parametersMap["by"];
        translation[2] = M_parametersMap["bz"];

        transformer.transformMesh(scaleVec, rotation, translation);

        R1 = computeRotationMatrix(0, rotation[0]);
        R2 = computeRotationMatrix(1, rotation[1]);
        R3 = computeRotationMatrix(2, rotation[2]);

        R = R1 * R2 * R3;
    }

    applyAffineTransformationGeometricFace(M_inlet,R,translation,scale);
    for (std::vector<GeometricFace>::iterator it = M_outlets.begin();
         it != M_outlets.end(); it++)
    {
        applyAffineTransformationGeometricFace(*it, R, translation,scale);
        it->print();
    }

    // Handle rotation along the axis of the inlet
    double angle = M_parametersMap["alpha_axis"];

    Matrix3D Raxis = computeRotationMatrix(M_inlet.M_normal,angle);
    Vector3D transZero = M_inlet.M_center - Raxis * M_inlet.M_center;

    auto foo = std::bind(rotationFunction,
                         std::placeholders::_1,
                         std::placeholders::_2,
                         std::placeholders::_3,
                         Raxis, transZero, 1.0);

    transformer.transformMesh(foo);

    for (std::vector<GeometricFace>::iterator it = M_outlets.begin();
         it != M_outlets.end(); it++)
    {
        applyAffineTransformationGeometricFace(*it, Raxis, transZero, 1.0);
    }
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
dumpMesh(std::string outdir, std::string meshdir, std::string outputName)
{
    if (!M_mesh)
    {
        std::string msg = "Mesh has not been read yet!\n";
        throw Exception(msg);
    }
    boost::filesystem::create_directory(outdir);

    GetPot exporterDatafile(meshdir + "datafiles/" + M_datafileName);
    LifeV::ExporterVTK<mesh_Type> exporter(exporterDatafile, outputName);
    exporter.setMeshProcId(M_mesh, M_comm->MyPID());

    FESpacePtr_Type dummyFespace(new FESpace_Type(M_mesh, "P1", 3, M_comm));
    vectorPtr_Type zero(new vector_Type(dummyFespace->map()) );
    zero->zero();

    exporter.addVariable(LifeV::ExporterData<mesh_Type>::ScalarField, "z",
                         dummyFespace, zero, 0);
    exporter.setPostDir(outdir);
    exporter.postProcess(0.0);
}

GeometricFace
BuildingBlock::
getOutlet(unsigned int outletIndex) const
{
    if (outletIndex >= M_outlets.size())
    {
        std::string msg = "Requesting acess to outlet that does not exist!";
        throw Exception(msg);
    }

    return M_outlets[outletIndex];
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
    M_inletRotationAxis = M_inletRotationAxis / M_inletRotationAxis.norm();
    M_inletAngle = std::acos(iNormal.dot(oNormal) /
                             (iNormal.norm() * oNormal.norm()));
}

void
BuildingBlock::
applyNonLinearTransformation()
{
}

void
BuildingBlock::
applyGlobalTransformation()
{
    applyAffineTransformation();
    applyNonLinearTransformation();
}

}  // namespace BuildingBlock
