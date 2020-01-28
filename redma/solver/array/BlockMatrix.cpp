#include "BlockMatrix.hpp"

namespace RedMA
{

template <>
template <>
void
BlockMatrix<MatrixEp>::
getRowProperty(std::shared_ptr<LifeV::MapEpetra>& outMap,
               const unsigned int& rowIndex) const
{
    outMap = nullptr;

    for (unsigned int j = 0; j < M_nCols; j++)
    {
        std::shared_ptr<LifeV::MapEpetra> curMap;
        block(rowIndex,j).getRowProperty(curMap);
        // we stop as soon as we find a non empty block
        if (curMap != nullptr)
        {
            outMap = curMap;
            return;
        }
    }
}

template <>
template <>
void
BlockMatrix<MatrixEp>::
getColProperty(std::shared_ptr<LifeV::MapEpetra>& outMap,
               const unsigned int& colIndex) const
{
    outMap = nullptr;

    for (unsigned int i = 0; i < M_nRows; i++)
    {
        std::shared_ptr<LifeV::MapEpetra> curMap;
        block(i,colIndex).getColProperty(curMap);
        // we stop as soon as we find a non empty block
        if (curMap != nullptr)
        {
            outMap = curMap;
            return;
        }
    }
}

template <>
template <>
void
BlockMatrix<MatrixEp>::
getRowsProperty(std::shared_ptr<LifeV::MapEpetra>& outMap) const
{
    outMap = nullptr;

    for (unsigned int i = 0; i < M_nRows; i++)
    {
        std::shared_ptr<LifeV::MapEpetra> curMap;
        getRowProperty(curMap, i);
        if (outMap == nullptr)
            outMap = curMap;
        else
            *outMap += *curMap;
    }
}

template <>
template <>
void
BlockMatrix<MatrixEp>::
getColsProperty(std::shared_ptr<LifeV::MapEpetra>& outMap) const
{
    outMap = nullptr;

    for (unsigned int j = 0; j < M_nCols; j++)
    {
        std::shared_ptr<LifeV::MapEpetra> curMap;
        getColProperty(curMap, j);
        if (outMap == nullptr)
            outMap = curMap;
        else
            *outMap += *curMap;
    }
}

template <>
template <>
void
BlockMatrix<MatrixEp>::
convertVectorType(const Epetra_MultiVector& inputVector,
                  BlockVector<VectorEp>& outputVector) const
{
    using namespace LifeV;

    std::shared_ptr<LifeV::MapEpetra> blockDomainMap;
    getColsProperty(blockDomainMap);
    VectorEpetra inputVectorEpetra(inputVector, blockDomainMap, Unique);

    outputVector.resize(M_nCols);

    unsigned int offset;
    for (unsigned int j = 0; j < M_nCols; j++)
    {
        std::shared_ptr<MapEpetra> curDomainMap;
        getColProperty(curDomainMap, j);

        std::shared_ptr<VectorEpetra> subVec(new VectorEpetra(curDomainMap, Unique));
        subVec->subset(inputVectorEpetra, *curDomainMap, offset, 0);

        outputVector.block(j) = subVec;

        offset += curDomainMap->mapSize();
    }

}

template <>
template <>
void
BlockMatrix<MatrixEp>::
convertVectorType(const BlockVector<VectorEp>& inputVector,
                  Epetra_MultiVector& outputVector) const
{
    using namespace LifeV;

    std::shared_ptr<LifeV::MapEpetra> blockDomainMap;
    getColsProperty(blockDomainMap);
    VectorEpetra outVectorEpetra(blockDomainMap, Unique);

    unsigned int offset;
    for (unsigned int j = 0; j < M_nCols; j++)
    {
        std::shared_ptr<MapEpetra> curDomainMap;
        getColProperty(curDomainMap, j);

        outVectorEpetra.subset(*inputVector.block(j).data(),
                               *curDomainMap, 0, offset);

        offset += curDomainMap->mapSize();
    }

    outputVector = dynamic_cast<Epetra_MultiVector&>(outVectorEpetra.epetraVector());

}

template <>
void
BlockMatrix<BlockMatrix<MatrixEp>>::
finalize()
{
    std::vector<unsigned int> rows(M_nRows);
    std::vector<unsigned int> cols(M_nCols);
    // resize empty blocks accordingly
    for (unsigned int i = 0; i < M_nRows; i++)
    {
        for (unsigned int j = 0; j < M_nCols; j++)
        {
            if (block(i,j).nRows() > 0)
            {
                rows[i] = block(i,j).nRows();
                cols[j] = block(i,j).nCols();
            }
        }
    }

    for (unsigned int i = 0; i < M_nRows; i++)
    {
        for (unsigned int j = 0; j < M_nCols; j++)
        {
            if (block(i,j).nRows() == 0)
            {
                block(i,j).resize(rows[i],rows[j]);
            }
        }
    }

    M_isFinalized = true;
}

}
