namespace RedMA
{

template <class InVectorType, class InMatrixType>
aAssembler<InVectorType, InMatrixType>::
aAssembler(const DataContainer& data) :
  M_data(data)
{
}

template <class InVectorType, class InMatrixType>
aAssembler<InVectorType, InMatrixType>::
aAssembler(const DataContainer& data, SHP(TreeNode) treeNode) :
  M_data(data),
  M_treeNode(treeNode)
{
}

}
