#include "LinearSystemSolver.hpp"

namespace RedMA
{

template <>
void
LinearSystemSolver<BlockVector<DenseVector>, BlockMatrix<DenseMatrix>>::
buildPreconditioner(const BlockMatrix<BlockMatrix<DenseMatrix>>& matrix)
{
}

template <>
void
LinearSystemSolver<BlockVector<DenseVector>, BlockMatrix<DenseMatrix>>::
solve(const BlockMatrix<BlockMatrix<DenseMatrix>>& matrix,
      const BlockVector<BlockVector<DenseVector>>& rhs,
      BlockVector<BlockVector<DenseVector>>& sol)
{
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

    M_numSolves++;
}

}
