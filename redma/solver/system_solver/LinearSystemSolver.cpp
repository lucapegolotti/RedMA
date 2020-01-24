#include "LinearSystemSolver.hpp"

namespace RedMA
{

template<>
BlockVector<VectorEp>
LinearSystemSolver<VectorEp, MatrixEp>::
solve(BlockMatrix<MatrixEp> matrix, BlockVector<VectorEp> rhs)
{
    BlockVector<VectorEp> sol;
    SHP(LinearOperatorEp) op;
    op->setup(matrix, nullptr);

    return sol;
}


}
