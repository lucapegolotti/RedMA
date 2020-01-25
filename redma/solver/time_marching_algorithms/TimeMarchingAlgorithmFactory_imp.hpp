namespace RedMA
{

template <class InVectorType, class InMatrixType>
SHP(aTimeMarchingAlgorithm<InVectorType COMMA InMatrixType>)
TimeMarchingAlgorithmFactory(const GetPot& datafile)
{
    std::shared_ptr<aTimeMarchingAlgorithm<InVectorType, InMatrixType>> ret;
    std::string algorithmString = datafile("time_discretization_algorithm","bdf");

    if (!std::strcmp(algorithmString.c_str(),"bdf"))
        ret.reset(new BDF<InVectorType, InMatrixType>(datafile));
    else
        throw new Exception("Requested time marching algorithm is not implemented!");

    return ret;
}

}
