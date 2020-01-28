#include "LinearSystemSolver.hpp"

namespace RedMA
{

// template<>
// BlockVector<VectorEp>
// LinearSystemSolver<VectorEp, MatrixEp>::
// solve(BlockMatrix<MatrixEp> matrix, BlockVector<VectorEp> rhs)
// {
//     BlockVector<VectorEp> sol;
//     SHP(LinearOperatorEp) op;
//     op->setup(matrix, nullptr);
//
//     return sol;
// }

template<>
BlockVector<BlockVector<VectorEp>>
LinearSystemSolver<BlockVector<VectorEp>, BlockMatrix<MatrixEp>>::
solve(BlockMatrix<BlockMatrix<MatrixEp>> matrix, BlockVector<BlockVector<VectorEp>> rhs)
{
    std::cout << "LinearSystemSolver" << std::endl << std::flush;
    BlockVector<BlockVector<VectorEp>> sol;
    LinearOperatorEp op;
    op.setup(matrix, nullptr);

    return sol;
}


}
