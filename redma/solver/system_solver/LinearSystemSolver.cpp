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

template<>
BlockVector<BlockVector<VectorEp>>
LinearSystemSolver<BlockVector<VectorEp>, BlockMatrix<MatrixEp>>::
solve(BlockMatrix<BlockMatrix<MatrixEp>> matrix, BlockVector<BlockVector<VectorEp>> rhs)
{
    std::cout << "solve1" << std::endl << std::flush;
    BlockVector<BlockVector<VectorEp>> sol;
    BlockMatrix<MatrixEp> matCollapsed;
    matrix.printPattern();
    std::cout << "----" << std::endl << std::flush;
    matrix.collapseBlocks(matCollapsed);
    matCollapsed.printPattern();
    std::cout << "----" << std::endl << std::flush;
    SHP(LinearOperatorEp) op;
    op->setup(matCollapsed, nullptr);
    std::cout << "solve2" << std::endl << std::flush;

    return sol;
}


}
