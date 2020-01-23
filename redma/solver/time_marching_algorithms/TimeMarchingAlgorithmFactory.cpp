#include "TimeMarchingAlgorithmFactory.hpp"

namespace RedMA
{

SHP(aTimeMarchingAlgorithm)
TimeMarchingAlgorithmFactory(const GetPot& datafile)
{
    SHP(aTimeMarchingAlgorithm) ret;
    std::string algorithmString = datafile("time_discretization_algorithm","bdf");

    if (!std::strcmp(algorithmString.c_str(),"bdf"))
        ret.reset(new BDF(datafile));
    else
        throw new Exception("Time Marching Algorithm is not implemented!");

    return ret;
}

}
