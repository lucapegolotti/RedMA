#include "LinearSystemSolver.hpp"

namespace RedMA
{

template<>
BlockVector<VectorEp>
LinearSystemSolver<VectorEp, MatrixEp>::
solve(BlockMatrix<MatrixEp> matrix, BlockVector<VectorEp> rh)
{
    BlockVector<VectorEp> sol;
    // magic of linear algebra
    LinearOperatorEp op;

    op.setup(matrix, nullptr);

    return sol;
}


}
