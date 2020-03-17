#include "SystemSolver.hpp"
namespace RedMA
{

template<class InVectorType, class InMatrixType>
SystemSolver<InVectorType, InMatrixType>::
SystemSolver(const DataContainer& data) :
  M_data(data),
  M_linearSystemSolver(data),
  M_isLinearProblem(false)
{
}

template<class InVectorType, class InMatrixType>
typename SystemSolver<InVectorType, InMatrixType>::BV
SystemSolver<InVectorType, InMatrixType>::
solve(FunctionFunctor<BV,BV> fun, FunctionFunctor<BV,BM> jac,
      BV initialGuess, int& status)
{
    BV sol = initialGuess;

    if (M_isLinearProblem)
    {
        BV incr(initialGuess.nRows());
        BV curFun = fun(sol);

        double err = curFun.norm2();

        std::ostringstream streamOb;
        streamOb << err;

        std::string msg = "[SystemSolver] linear solve,";
        msg += " residual = " + streamOb.str() + "\n";
        printlog(GREEN, msg, M_data.getVerbose());

        incr.zero();
        BM curJac = jac(sol);
        curJac.finalize();

        M_linearSystemSolver.solve(curJac, curFun, incr);

        sol -= incr;

        status = 0;
    }
    else
    {
        M_solverStatistics.resize(0);

        BV incr(initialGuess.nRows());
        BV curFun;
        BV sol = initialGuess;

        double tol = M_data("newton_method/tol", 1e-5);
        unsigned int maxit = M_data("newton_method/maxit", 10);

        std::string msg;

        double err = tol + 1;
        unsigned int count = 0;
        double initialError = 1;
        while (err / initialError > tol && count < maxit)
        {
            if (count == 0)
            {
                curFun = fun(sol);
                err = curFun.norm2();
                initialError = err;
            }

            std::ostringstream streamOb;
            streamOb << err / initialError;

            msg = "[SystemSolver] Newtons algorithm,";
            msg += " rel error = " + streamOb.str();
            streamOb.str("");
            streamOb.clear();
            streamOb << err;
            msg += ", abs error = " + streamOb.str();
            msg += ", iteration = " + std::to_string(count+1) + "\n";
            printlog(GREEN, msg, M_data.getVerbose());

            if (err / initialError > tol)
            {
                incr.zero();
                BM curJac = jac(sol);
                curJac.finalize();
                M_linearSystemSolver.solve(curJac, curFun, incr);
                M_solverStatistics.push_back(M_linearSystemSolver.getSolverStatistics());
            }
            sol -= incr;
            count++;

            curFun = fun(sol);
            err = curFun.norm2();
        }

        // newton method has failed
        if (count == maxit && err / initialError > tol)
        {
            status = -1;

            std::ostringstream streamOb;

            streamOb << err / initialError;

            msg = "[SystemSolver] Newtons algorithm, failed";
            msg += " rel error = " + streamOb.str();
            streamOb.str("");
            streamOb.clear();
            streamOb << err;
            msg += ", abs error = " + streamOb.str();
            msg += ", iteration = " + std::to_string(count+1) + "\n";
            printlog(GREEN, msg, M_data.getVerbose());
        }
        else
        {
            status = 0;

            std::ostringstream streamOb;

            streamOb << err / initialError;

            msg = "[SystemSolver] Newtons algorithm converged,";
            msg += " rel error = " + streamOb.str();
            streamOb.str("");
            streamOb.clear();
            streamOb << err;
            msg += ", abs error = " + streamOb.str();
            msg += ", iteration = " + std::to_string(count+1) + "\n";

            printlog(GREEN, msg, M_data.getVerbose());
        }
    }
    return sol;
}

}
