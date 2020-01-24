namespace RedMA
{

template<class InVectorType, class InMatrixType>
LinearSystemSolver<InVectorType, InMatrixType>::
LinearSystemSolver(const GetPot& datafile) :
  M_datafile(datafile)
{
}

template<class InVectorType, class InMatrixType>
typename LinearSystemSolver<InVectorType, InMatrixType>::BV
LinearSystemSolver<InVectorType, InMatrixType>::
solve(BM matrix, BV rh)
{
    BV sol;
    // magic of linear algebra
    return sol;
}

}
