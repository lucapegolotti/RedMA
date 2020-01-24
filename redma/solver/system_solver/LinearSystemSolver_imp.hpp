namespace RedMA
{

template<class InVectorType, class InMatrixType>
LinearSystemSolver<InVectorType, InMatrixType>::
LinearSystemSolver(const GetPot& datafile) :
  M_datafile(datafile)
{
}

}
