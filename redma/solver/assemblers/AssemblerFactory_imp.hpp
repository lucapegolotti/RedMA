namespace RedMA
{

template<class InVectorType, class InMatrixType>
SHP(aAssembler<InVectorType COMMA InMatrixType>)
AssemblerFactory(const DataContainer& data, SHP(TreeNode) treeNode)
{
    std::shared_ptr<aAssembler<InVectorType, InMatrixType>> ret;
    std::string assemblerString = data("assembler/type","stokes");

    if (!std::strcmp(assemblerString.c_str(),"stokes"))
        ret.reset(new StokesAssembler<InVectorType,InMatrixType>(data,
                                                                 treeNode));
    else
        throw new Exception("Time Marching Algorithm is not implemented!");

    return ret;
}

}
