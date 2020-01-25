namespace RedMA
{

template <class InVectorType, class InMatrixType>
BlockAssembler<InVectorType, InMatrixType>::
BlockAssembler(const GetPot& datafile, SHP(TreeStructure) tree) :
  aAssembler<InVectorType, InMatrixType>(datafile),
  M_tree(tree)
{

}

template <class InVectorType, class InMatrixType>
void
BlockAssembler<InVectorType, InMatrixType>::
exportSolution(const double& t)
{

}

template <class InVectorType, class InMatrixType>
void
BlockAssembler<InVectorType, InMatrixType>::
postProcess()
{

}

template <class InVectorType, class InMatrixType>
BlockMatrix<InMatrixType>
BlockAssembler<InVectorType, InMatrixType>::
getMass(const double& time, const BlockVector<InVectorType>& sol)
{

}

template <class InVectorType, class InMatrixType>
BlockVector<InVectorType>
BlockAssembler<InVectorType, InMatrixType>::
getRightHandSide(const double& time, const BlockVector<InVectorType>& sol)
{

}

template <class InVectorType, class InMatrixType>
BlockMatrix<InMatrixType>
BlockAssembler<InVectorType, InMatrixType>::
getJacobianRightHandSide(const double& time, const BlockVector<InVectorType>& sol)
{

}

}
