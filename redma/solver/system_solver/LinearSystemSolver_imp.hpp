namespace RedMA
{

template<class InVectorType, class InMatrixType>
LinearSystemSolver<InVectorType, InMatrixType>::
LinearSystemSolver(const DataContainer& data) :
  M_data(data),
  M_numSolves(0)
{
}

}
