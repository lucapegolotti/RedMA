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
        unsigned int recomputeevery = M_data("preconditioner/recomputeevery", 1);
        if (M_prec == nullptr || (M_numSolves % recomputeevery) == 0)
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

    Chrono chrono;
    chrono.start();
    printlog(MAGENTA, "[LinearSystemSolver] solve ...", M_data.getVerbose());
    M_statistics.M_numIterations = M_invOper->invert(rhs, sol);

    M_statistics.M_solveTime = chrono.diff();
    std::string msg = "done, in ";
    msg += std::to_string(M_statistics.M_solveTime);
    msg += " seconds\n";
    printlog(GREEN, msg, M_data.getVerbose());

    M_numSolves++;
}


template <>
void
LinearSystemSolver<VectorEp, MatrixEp>::
solve(const BlockMatrix<MatrixEp>& matrix,
      const BlockVector<VectorEp>& rhs,
      BlockVector<VectorEp>& sol)
{
    throw new Exception("This specialization of LinearSystemSolver::solve "
                        "should not be used");
}

template <>
void
LinearSystemSolver<VectorEp, MatrixEp>::
computeSchurComplementDense(const BlockMatrix<MatrixEp>& matrix)
{
    throw new Exception("This specialization of LinearSystemSolver::computeSchurComplementDense "
                        "should not be used");
}



}
