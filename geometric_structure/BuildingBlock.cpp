#include "BuildingBlock.hpp"
#include <iostream>

namespace ReMA
{

BuildingBlock::BuildingBlock(commPtr_Type comm, bool verbose) :
  M_comm(comm),
  M_verbose(verbose)
{
    if (M_comm->MyPID() != 0)
        M_verbose = false;
}

void BuildingBlock::setParameterValue(std::string key, double value)
{
    if (M_parametersMap.find(key) != M_parametersMap.end())
    {
        M_parametersMap[key] = value;
    }
    else
    {
        std::string errorMsg = "Parameter with key " + key + " and value " +
                                std::to_string(value) + " is not contained in " +
                                M_name + " building block!\n";

        throw Exception(errorMsg);
    }
}

int BuildingBlock::readMesh(std::string meshdir)
{
    printlog(GREEN, "[BuildingBlock] reading mesh ...\n", M_verbose);

    meshPtr_Type fullMesh(new mesh_Type(M_comm));
    LifeV::MeshData meshData;
    GetPot meshDatafile(meshdir + "datafiles/" + M_datafileName);
    meshData.setup(meshDatafile, "mesh");
    meshData.setMeshDir(meshdir);

    LifeV::readMesh(*fullMesh,meshData);

    LifeV::MeshPartitioner<mesh_Type> meshPart;

    // small trick to redirect std cout
    std::streambuf* curBuf = std::cout.rdbuf();
    std::ostringstream strCout;
    std::cout.rdbuf(strCout.rdbuf());
    meshPart.doPartition(fullMesh, M_comm);
    std::cout.rdbuf(curBuf);
    printlog(YELLOW, strCout.str(), M_verbose);

    printlog(GREEN, "done\n", M_verbose);

    return 0;
}

}  // namespace BuildingBlock
