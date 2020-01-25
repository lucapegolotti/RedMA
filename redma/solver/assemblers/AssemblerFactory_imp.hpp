namespace RedMA
{

template<class InVectorType, class InMatrixType>
SHP(aAssembler<InVectorType COMMA InMatrixType>)
AssemblerFactory(const GetPot& datafile, SHP(BuildingBlock) buildingBlock)
{
    std::shared_ptr<aAssembler<InVectorType, InMatrixType>> ret;
    std::string assemblerString = datafile("assembler/type","stokes");

    if (!std::strcmp(assemblerString.c_str(),"stokes"))
        ret.reset(new StokesAssembler<InVectorType,InMatrixType>(datafile,
                                                                 buildingBlock));
    else
        throw new Exception("Time Marching Algorithm is not implemented!");

    return ret;
}

}
