#include "SystemSolver.hpp"
namespace RedMA
{

template<class InVectorType, class InMatrixType>
SystemSolver<InVectorType, InMatrixType>::
SystemSolver(const DataContainer& data) :
  M_data(data),
  M_linearSystemSolver(data)
{
}

template<class InVectorType, class InMatrixType>
typename SystemSolver<InVectorType, InMatrixType>::BV
SystemSolver<InVectorType, InMatrixType>::
solve(FunctionFunctor<BV,BV> fun, FunctionFunctor<BV,BM> jac,
      BV initialGuess, int& status)
{
    BV incr;
    BV curFun;
    BV sol = initialGuess;

    double tol = M_data("newton_method/tol", 1e-5);
    unsigned int maxit = M_data("newton_method/maxit", 10);

    double err = tol + 1;
    unsigned int count = 0;
    while (err > tol && count < maxit)
    {
        curFun = fun(sol);
        err = curFun.norm2();

        if (err > tol)
        {
            incr.zero();
            BM curJac = jac(sol);

            incr = M_linearSystemSolver.solve(curJac, curFun);
        }
        sol -= incr;
        count++;
    }
    // newton method has failed
    if (count == maxit)
        status = -1;

    status = 0;
    return sol;
}

}
