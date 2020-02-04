#include "LinearSystemSolver.hpp"

namespace RedMA
{

template <>
void
LinearSystemSolver<BlockVector<VectorEp>, BlockMatrix<MatrixEp>>::
buildPreconditioner(const BlockMatrix<BlockMatrix<MatrixEp>>& matrix)
{
    std::string precType = M_data("preconditioner/outer", "saddlepoint");

    if (!std::strcmp(precType.c_str(), "saddlepoint"))
    {
        M_prec.reset(new SaddlePointPreconditionerEp(M_data, matrix));
    }
    else
    {
        throw new Exception("Unkown type of preconditioner!");
    }
    M_statistics.M_precSetupTime = M_prec->getSetupTime();
}

template <>
void
LinearSystemSolver<BlockVector<VectorEp>, BlockMatrix<MatrixEp>>::
solve(const BlockMatrix<BlockMatrix<MatrixEp>>& matrix,
      const BlockVector<BlockVector<VectorEp>>& rhs,
      BlockVector<BlockVector<VectorEp>>& sol)
{
    M_oper.reset(new LinearOperatorEp(matrix));

    M_invOper.reset(new InverseOperatorEp(M_data));
    M_invOper->setOperator(M_oper);

    buildPreconditioner(matrix);

    M_invOper->setPreconditioner(M_prec);

    LifeV::LifeChrono chrono;
    chrono.start();
    printlog(MAGENTA, "[LinearSystemSolver] solve ...", M_data.getVerbose());

    M_statistics.M_numIterations = M_invOper->invert(rhs, sol);

    M_statistics.M_solveTime = chrono.diff();
    std::string msg = "done, in ";
    msg += std::to_string(M_statistics.M_solveTime);
    msg += " seconds\n";
    printlog(GREEN, msg, M_data.getVerbose());
}

template <>
void
LinearSystemSolver<Double, Double>::
solve(const BlockMatrix<Double>& matrix,
      const BlockVector<Double>& rhs,
      BlockVector<Double>& sol)
{
    if (matrix.nRows() > 1 || matrix.nCols() > 1 || rhs.nRows() > 1)
        throw new Exception("solve with blocks of doubles is implemented only in"
                            " dimension 1");

    sol.block(0).data() = rhs.block(0).data() / matrix.block(0,0).data();
}


}
