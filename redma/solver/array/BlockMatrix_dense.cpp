#include "BlockMatrix.hpp"

namespace RedMA
{

template <>
template <>
void
BlockMatrix<DenseMatrix>::
getRowProperty(std::shared_ptr<LifeV::MapEpetra>& outMap,
               const unsigned int& rowIndex) const
{
}

template <>
template <>
void
BlockMatrix<DenseMatrix>::
getColProperty(std::shared_ptr<LifeV::MapEpetra>& outMap,
               const unsigned int& colIndex) const
{
}

template <>
template <>
void
BlockMatrix<DenseMatrix>::
getRowsProperty(std::shared_ptr<LifeV::MapEpetra>& outMap) const
{
}

template <>
template <>
void
BlockMatrix<DenseMatrix>::
getColsProperty(std::shared_ptr<LifeV::MapEpetra>& outMap) const
{
}

template <>
template <>
void
BlockMatrix<DenseMatrix>::
convertVectorType(const Epetra_MultiVector& inputVector,
                  BlockVector<VectorEp>& outputVector) const
{
}

template <>
template <>
void
BlockMatrix<DenseMatrix>::
convertVectorType(const BlockVector<VectorEp>& inputVector,
                  Epetra_MultiVector& outputVector) const
{
}

template <>
void
BlockMatrix<DenseMatrix>::
finalize(BlockMatrix<BlockMatrix<DenseMatrix>>* father,
         unsigned int* myRow, unsigned int* myCol)
{
}

template <>
void
BlockMatrix<BlockMatrix<DenseMatrix>>::
finalize(BlockMatrix<BlockMatrix<BlockMatrix<DenseMatrix>>>* father,
         unsigned int* myRow, unsigned int* myCol)
{
}

template <>
void
BlockMatrix<DenseMatrix>::
printPattern() const
{
    for (unsigned int i = 0; i < M_nRows; i++)
    {
        for (unsigned int j = 0; j < M_nCols; j++)
        {
            auto curblock = block(i,j);
            if (curblock.data())
                std::cout << "(" << curblock.data()->M() <<
                             "," << curblock.data()->N() << ")";
            else
                std::cout << "o";
            std::cout << "\t\t";
        }
        std::cout << "\n";
    }
}

}
