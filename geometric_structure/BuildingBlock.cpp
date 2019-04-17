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
  M_verbose(verbose)
{
    if (M_comm->MyPID() != 0)
        M_verbose = false;

    // rotation angle
    M_parametersMap["alphax"] = 0.0;
    M_parametersMap["alphay"] = 0.0;
    M_parametersMap["alphaz"] = 0.0;

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

void
BuildingBlock::
applyAffineTransformation()
{
    LifeV::MeshUtility::MeshTransformer<mesh_Type> transformer(*M_mesh);

    Vector3D scale(M_parametersMap["scale"],
                   M_parametersMap["scale"],
                   M_parametersMap["scale"]);

    Vector3D rotation(-M_parametersMap["alphax"],
                      -M_parametersMap["alphay"],
                      -M_parametersMap["alphaz"]);

    Vector3D translation(M_parametersMap["bx"],
                         M_parametersMap["by"],
                         M_parametersMap["bz"]);

    transformer.transformMesh(scale, rotation, translation);

    Matrix3D R, R1, R2, R3, S;

    R1 = computeRotationMatrix(0, rotation[0]);
    R2 = computeRotationMatrix(1, rotation[1]);
    R3 = computeRotationMatrix(2, rotation[2]);

    S(0,0) = scale[0];
    S(0,1) = 0.;
    S(0,2) = 0.;
    S(1,0) = 0.;
    S(1,1) = scale[1];
    S(1,2) = 0.;
    S(2,0) = 0.;
    S(2,1) = 0.;
    S(2,2) = scale[2];

    R = R3 * S;
    R = R2 * R;
    R = R1 * R;

    applyAffineTransformationGeometricFace(M_inlet,R,translation,scale[0]);
    for (std::vector<GeometricFace>::iterator it = M_outlets.begin();
         it != M_outlets.end(); it++)
    {
        applyAffineTransformationGeometricFace(*it,R,translation,scale[0]);
    }
    std::cout << "Apply affine transf\n" << std::endl;
    M_inlet.print();
}

void
BuildingBlock::
applyAffineTransformationGeometricFace(GeometricFace& face,
                                       const Matrix3D& affineMatrix,
                                       const Vector3D& translation,
                                       const double& scale)
{
    face.M_center = affineMatrix * face.M_center;
    face.M_center = face.M_center + translation;

    face.M_normal = affineMatrix * face.M_normal;

    face.M_radius = face.M_radius * scale;
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
    parentOutlet.print();

    M_parametersMap["bx"] = parentOutlet.M_center[0] - M_inlet.M_center[0];
    M_parametersMap["by"] = parentOutlet.M_center[1] - M_inlet.M_center[1];
    M_parametersMap["bz"] = parentOutlet.M_center[2] - M_inlet.M_center[2];

    M_parametersMap["scale"] = parentOutlet.M_radius / M_inlet.M_radius;

    unsigned int coors1[3] = {1,2,0};
    unsigned int coors2[3] = {2,0,1};

    std::vector<double> angles;
    std::vector<double> dets;

    Vector3D iNormal = -1 * M_inlet.M_normal;
    Vector3D oNormal = parentOutlet.M_normal;

    // here we compute each angle by considering the projection of the vectors
    // on the remaining plane (e.g. alphax->projection on yz)
    for (int i = 0; i < 3; i++)
    {

        unsigned int i1 = coors1[i];
        unsigned int i2 = coors2[i];

        double norminl = std::sqrt(iNormal[i1] *
                                   iNormal[i1] +
                                   iNormal[i2] *
                                   iNormal[i2]);

        double normout = std::sqrt(oNormal[i1] *
                                   oNormal[i1] +
                                   oNormal[i2] *
                                   oNormal[i2]);

        double val;
        if (norminl < 1e-14 || normout < 1e-14)
        {
            val = 1;
            dets.push_back(0.0);
        }
        else
        {
            val = (iNormal[i1] * oNormal[i1] +
                   iNormal[i2] * oNormal[i2]) /
                   (norminl * normout);
            double determ = iNormal[i1] * oNormal[i2] -
                            iNormal[i2] * oNormal[i1];

            dets.push_back(determ);
        }
        double angle;
        if (dets[i] == 0)
            angle = 0;
        else if (dets[i] > 0)
            angle = -std::acos(val) + M_PI;
        else if (dets[i] < 0)
            angle = std::acos(val);

        Matrix3D R = computeRotationMatrix(i, angle);
        oNormal = R * oNormal;

        angles.push_back(angle);
    }

    M_parametersMap["alphax"] = angles[0];
    M_parametersMap["alphay"] = angles[1];
    M_parametersMap["alphaz"] = angles[2];
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
