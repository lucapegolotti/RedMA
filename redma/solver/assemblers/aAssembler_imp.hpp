namespace RedMA
{

template <class InVectorType, class InMatrixType>
aAssembler<InVectorType, InMatrixType>::
aAssembler(const GetPot& datafile) :
  M_datafile(datafile)
{
}

template <class InVectorType, class InMatrixType>
aAssembler<InVectorType, InMatrixType>::
aAssembler(const GetPot& datafile, SHP(TreeNode) treeNode) :
  M_datafile(datafile),
  M_treeNode(treeNode)
{
}

}
