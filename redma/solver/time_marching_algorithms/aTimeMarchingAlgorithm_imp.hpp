namespace RedMA
{

template<class InVectorType, class InMatrixType>
aTimeMarchingAlgorithm<InVectorType,InMatrixType>::
aTimeMarchingAlgorithm(const DataContainer& data,
                       SHP(aAssembler<InVectorType COMMA InMatrixType>) assembler) :
  M_data(data),
  M_systemSolver(data),
  M_assembler(assembler)
{
}

}
