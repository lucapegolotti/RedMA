#include "TimeMarchingAlgorithmFactory.hpp"

namespace RedMA
{

SHP(aTimeMarchingAlgorithm)
TimeMarchingAlgorithmFactory(const DataContainer& data)
{
    SHP(aTimeMarchingAlgorithm) ret;
    std::string algorithmString = data("time_discretization/algorithm","bdf");

    if (!std::strcmp(algorithmString.c_str(),"bdf"))
        ret.reset(new BDF(data));
    else if (!std::strcmp(algorithmString.c_str(),"alpha"))
        ret.reset(new GeneralizedAlphaMethod(data));
    // else if (!std::strcmp(algorithmString.c_str(),"alpha1storderp"))
    //     ret.reset(new GeneralizedAlphaMethod<InVectorType, InMatrixType>(data));
    else
        throw new Exception("Requested time marching algorithm is not implemented!");

    return ret;
}

SHP(aTimeMarchingAlgorithm)
TimeMarchingAlgorithmFactory(const DataContainer& data,
                             SHP(aFunctionProvider) funProvider)
{
    SHP(aTimeMarchingAlgorithm) ret;
    std::string algorithmString = data("time_discretization/algorithm","bdf");

    if (!std::strcmp(algorithmString.c_str(),"bdf"))
        ret.reset(new BDF(data, funProvider));
    else if (!std::strcmp(algorithmString.c_str(),"alpha"))
        ret.reset(new GeneralizedAlphaMethod(data, funProvider));
    // else if (!std::strcmp(algorithmString.c_str(),"alpha1storderp"))
    //     ret.reset(new GeneralizedAlphaMethod1stOrderPressure<InVectorType, InMatrixType>(data, funProvider));
    else
        throw new Exception("Requested time marching algorithm is not implemented!");

    return ret;
}

}
