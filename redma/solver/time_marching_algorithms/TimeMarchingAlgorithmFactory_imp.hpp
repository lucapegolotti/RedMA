namespace RedMA
{

template <class InVectorType, class InMatrixType>
SHP(aTimeMarchingAlgorithm<InVectorType COMMA InMatrixType>)
TimeMarchingAlgorithmFactory(const DataContainer& data)
{
    std::shared_ptr<aTimeMarchingAlgorithm<InVectorType, InMatrixType>> ret;
    std::string algorithmString = data("time_discretization/algorithm","bdf");

    if (!std::strcmp(algorithmString.c_str(),"bdf"))
        ret.reset(new BDF<InVectorType, InMatrixType>(data));
    else if (!std::strcmp(algorithmString.c_str(),"alpha"))
        ret.reset(new GeneralizedAlphaMethod<InVectorType, InMatrixType>(data));
    // else if (!std::strcmp(algorithmString.c_str(),"alpha1storderp"))
    //     ret.reset(new GeneralizedAlphaMethod<InVectorType, InMatrixType>(data));
    else
        throw new Exception("Requested time marching algorithm is not implemented!");

    return ret;
}

template <class InVectorType, class InMatrixType>
SHP(aTimeMarchingAlgorithm<InVectorType COMMA InMatrixType>)
TimeMarchingAlgorithmFactory(const DataContainer& data,
                             SHP(aFunctionProvider<InVectorType COMMA InMatrixType>) funProvider)
{
    std::shared_ptr<aTimeMarchingAlgorithm<InVectorType, InMatrixType>> ret;
    std::string algorithmString = data("time_discretization/algorithm","bdf");

    if (!std::strcmp(algorithmString.c_str(),"bdf"))
        ret.reset(new BDF<InVectorType, InMatrixType>(data, funProvider));
    else if (!std::strcmp(algorithmString.c_str(),"alpha"))
        ret.reset(new GeneralizedAlphaMethod<InVectorType, InMatrixType>(data, funProvider));
    // else if (!std::strcmp(algorithmString.c_str(),"alpha1storderp"))
    //     ret.reset(new GeneralizedAlphaMethod1stOrderPressure<InVectorType, InMatrixType>(data, funProvider));
    else
        throw new Exception("Requested time marching algorithm is not implemented!");

    return ret;
}

}
