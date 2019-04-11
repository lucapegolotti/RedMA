#include "BuildingBlock.hpp"
#include <iostream>

namespace ReMA
{

BuildingBlock::BuildingBlock()
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

}  // namespace BuildingBlock
