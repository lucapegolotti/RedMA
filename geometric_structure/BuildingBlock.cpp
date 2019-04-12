#include "BuildingBlock.hpp"
#include <iostream>

namespace ReMA
{

BuildingBlock::BuildingBlock(commPtr_Type comm, bool verbose) :
  M_comm(comm),
  M_verbose(verbose)
{
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

void BuildingBlock::readMesh()
{
    if (M_verbose)
        printlog(GREEN,"[BuildingBlock] reading mesh ...");

    meshPtr_Type fullMesh(new mesh_Type(M_comm));

}

}  // namespace BuildingBlock
