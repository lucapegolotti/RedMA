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

    M_invOper->invert(rhs, sol);
}


}
