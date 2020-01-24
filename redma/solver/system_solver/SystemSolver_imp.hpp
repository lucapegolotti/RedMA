#include "SystemSolver.hpp"
namespace RedMA
{

template<class InVectorType, class InMatrixType>
typename SystemSolver<InVectorType, InMatrixType>::BV
SystemSolver<InVectorType, InMatrixType>::
solve(FunctionFunctor<BV,BV> function,
      FunctionFunctor<BV,BM> jacobian)
{

}

}
