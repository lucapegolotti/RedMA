#include "TimeMarchingAlgorithmFactory.hpp"

namespace RedMA
{

shp<aTimeMarchingAlgorithm>
TimeMarchingAlgorithmFactory(const DataContainer& data)
{
    shp<aTimeMarchingAlgorithm> ret;
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

shp<aTimeMarchingAlgorithm>
TimeMarchingAlgorithmFactory(const DataContainer& data,
                             shp<aFunctionProvider> funProvider)
{
    shp<aTimeMarchingAlgorithm> ret;
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
