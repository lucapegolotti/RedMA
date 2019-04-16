#include "BuildingBlock.hpp"
#include <iostream>

namespace RedMA
{

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
}

}  // namespace BuildingBlock
