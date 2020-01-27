namespace RedMA
{

template<class InVectorType, class InMatrixType>
aTimeMarchingAlgorithm<InVectorType,InMatrixType>::
aTimeMarchingAlgorithm(const DataContainer& data) :
  M_data(data),
  M_systemSolver(data)
{
}

}
