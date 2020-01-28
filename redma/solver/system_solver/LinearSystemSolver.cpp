#include "LinearSystemSolver.hpp"

namespace RedMA
{

template<>
void
LinearSystemSolver<BlockVector<VectorEp>, BlockMatrix<MatrixEp>>::
solve(const BlockMatrix<BlockMatrix<MatrixEp>>& matrix,
      const BlockVector<BlockVector<VectorEp>>& rhs,
      BlockVector<BlockVector<VectorEp>>& sol)
{
    M_oper.reset(new LinearOperatorEp(matrix));

    M_invOper.reset(new InverseOperatorEp(M_data));
    M_invOper->setOperator(M_oper);
    // M_invOper->setPreconditioner(nullptr);

    M_invOper->invert(rhs, sol);
}

}
