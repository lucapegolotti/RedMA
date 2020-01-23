namespace RedMA
{

template <class DataType>
SHP(aTimeMarchingAlgorithm<DataType>)
TimeMarchingAlgorithmFactory(const GetPot& datafile)
{
    SHP(aTimeMarchingAlgorithm<DataType>) ret;
    std::string algorithmString = datafile("time_discretization_algorithm","bdf");

    if (!std::strcmp(algorithmString.c_str(),"bdf"))
        ret.reset(new BDF<DataType>(datafile));
    else
        throw new Exception("Time Marching Algorithm is not implemented!");

    return ret;
}

}
