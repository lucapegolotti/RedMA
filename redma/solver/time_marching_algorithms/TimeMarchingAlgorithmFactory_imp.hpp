namespace RedMA
{

template <class InVectorType, class InMatrixType>
SHP(aTimeMarchingAlgorithm<InVectorType COMMA InMatrixType>)
TimeMarchingAlgorithmFactory(const DataContainer& data,
                             SHP(aAssembler<InVectorType COMMA InMatrixType>) assembler)
{
    std::shared_ptr<aTimeMarchingAlgorithm<InVectorType, InMatrixType>> ret;
    std::string algorithmString = data("time_discretization_algorithm","bdf");

    if (!std::strcmp(algorithmString.c_str(),"bdf"))
        ret.reset(new BDF<InVectorType, InMatrixType>(data, assembler));
    else
        throw new Exception("Requested time marching algorithm is not implemented!");

    return ret;
}

}
