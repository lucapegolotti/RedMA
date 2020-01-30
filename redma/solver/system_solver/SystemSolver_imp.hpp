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

    std::string msg;

    double err = tol + 1;
    unsigned int count = 0;
    while (err > tol && count < maxit)
    {
        if (count == 0)
        {
            curFun = fun(sol);
            err = curFun.norm2();
        }

        std::ostringstream streamOb;
        streamOb << err;

        msg = "[SystemSolver] Newtons algorithm,";
        msg += " error = " + streamOb.str();
        msg += ", iteration = " + std::to_string(count+1) + "\n";
        printlog(GREEN, msg, M_data.getVerbose());

        if (err > tol)
        {
            incr.zero();
            BM curJac = jac(sol);
            curJac.finalize();
            M_linearSystemSolver.solve(curJac, curFun, incr);
        }
        sol -= incr;
        count++;

        curFun = fun(sol);
        err = curFun.norm2();
    }

    std::cout << "err == " << err << std::endl << std::flush;
    // newton method has failed
    if (count == maxit && err > tol)
        status = -1;
    else
    {
        std::ostringstream streamOb;
        streamOb << err;

        msg = "[SystemSolver] Newtons algorithm, convergence,";
        msg += " error = " + streamOb.str();
        msg += ", iteration = " + std::to_string(count+1) + "\n";
        printlog(GREEN, msg, M_data.getVerbose());
    }

    status = 0;
    return sol;
}

}
