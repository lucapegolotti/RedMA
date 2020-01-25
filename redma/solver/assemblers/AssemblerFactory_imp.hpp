namespace RedMA
{

template<class InVectorType, class InMatrixType>
SHP(aAssembler<InVectorType COMMA InMatrixType>)
AssemblerFactory(const GetPot& datafile, SHP(TreeNode) treeNode)
{
    std::shared_ptr<aAssembler<InVectorType, InMatrixType>> ret;
    std::string assemblerString = datafile("assembler/type","stokes");

    if (!std::strcmp(assemblerString.c_str(),"stokes"))
        ret.reset(new StokesAssembler<InVectorType,InMatrixType>(datafile,
                                                                 treeNode));
    else
        throw new Exception("Time Marching Algorithm is not implemented!");

    return ret;
}

}
