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

void
BuildingBlock::
applyAffineTransformation()
{
    LifeV::MeshUtility::MeshTransformer<mesh_Type> transformer(*M_mesh);

    Vector3D scale(M_parametersMap["scale"],
                   M_parametersMap["scale"],
                   M_parametersMap["scale"]);

    Vector3D rotation(M_parametersMap["alphax"],
                      M_parametersMap["alphay"],
                      M_parametersMap["alphaz"]);

    Vector3D translation(M_parametersMap["bx"],
                         M_parametersMap["by"],
                         M_parametersMap["bz"]);

    transformer.transformMesh(scale, rotation, translation);

    Matrix3D R, R1, R2, R3, S;

    R1(0,0) = 1.;
    R1(0,1) = 0.;
    R1(0,2) = 0.;
    R1(1,0) = 0.;
    R1(1,1) = std::cos (rotation[0]);
    R1(1,2) = -std::sin (rotation[0]);
    R1(2,0) = 0.;
    R1(2,1) = std::sin (rotation[0]);
    R1(2,2) = std::cos (rotation[0]);

    R2(0,0) = std::cos (rotation[1]);
    R2(0,1) = 0.;
    R2(0,2) = std::sin (rotation[1]);
    R2(1,0) = 0.;
    R2(1,1) = 1.;
    R2(1,2) = 0.;
    R2(2,0) = -std::sin (rotation[1]);
    R2(2,1) = 0.;
    R2(2,2) = std::cos (rotation[1]);

    R3(0,0) = std::cos (rotation[2]);
    R3(0,1) = -std::sin (rotation[2]);
    R3(0,2) = 0.;
    R3(1,0) = std::sin (rotation[2]);
    R3(1,1) = std::cos (rotation[2]);
    R3(1,2) = 0.;
    R3(2,0) = 0.;
    R3(2,1) = 0.;
    R3(2,2) = 1.;

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
    GetPot exporterDatafile(meshdir + "datafiles/" + M_datafileName);
    LifeV::ExporterVTK<mesh_Type> exporter(exporterDatafile, outputName);
    exporter.setMeshProcId(M_mesh, M_comm->MyPID());

    FESpacePtr_Type dummyFespace(new FESpace_Type(M_mesh, "P1", 3, M_comm));
    vectorPtr_Type zero( new vector_Type(dummyFespace->map()) );
    zero->zero();

    exporter.addVariable(LifeV::ExporterData<mesh_Type>::ScalarField, "z",
                         dummyFespace, zero, 0);
    exporter.postProcess(0.0);
}

}  // namespace BuildingBlock
