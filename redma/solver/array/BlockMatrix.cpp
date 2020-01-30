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
BlockMatrix<MatrixEp>::
finalize(BlockMatrix<BlockMatrix<MatrixEp>>* father,
         unsigned int* myRow, unsigned int* myCol)
{
    typedef LifeV::MapEpetra                MP;
    typedef std::shared_ptr<MP>             MPptr;
    std::vector<MPptr> ranges(M_nRows);
    std::vector<MPptr> domains(M_nCols);

    if (!father || !myRow || !myCol)
    {
        for (unsigned int i = 0; i < M_nRows; i++)
        {
            for (unsigned int j = 0; j < M_nCols; j++)
            {
                if (block(i,j).data())
                {
                    ranges[i].reset(new MP(block(i,j).data()->rangeMap()));
                    domains[j].reset(new MP(block(i,j).data()->domainMap()));
                }
            }
        }
    }
    else
    {
        for (unsigned int jouter = 0; jouter < father->nCols(); jouter++)
        {
            auto curblock = father->block(*myRow,jouter);

            if (curblock.nRows() > 0)
            {
                if (curblock.nRows() != M_nRows)
                    throw new Exception("Error in finalize. Dimensions not consistent");

                for (unsigned int i = 0; i <curblock.nRows(); i++)
                {
                    for (unsigned int j = 0; j < curblock.nCols(); j++)
                    {
                        if (curblock.block(i,j).data())
                        {
                            ranges[i].reset(new MP(curblock.block(i,j).data()->rangeMap()));
                            break;
                        }
                    }
                }
            }
        }

        for (unsigned int iouter = 0; iouter < father->nRows(); iouter++)
        {
            auto curblock = father->block(iouter,*myCol);
            if (curblock.nRows() > 0)
            {
                if (curblock.nCols() != M_nCols)
                    throw new Exception("Error in finalize. Dimensions not consistent");

                for (unsigned int j = 0; j < curblock.nCols(); j++)
                {
                    for (unsigned int i = 0; i < curblock.nRows(); i++)
                    {
                        if (curblock.block(i,j).data())
                        {
                            domains[j].reset(new MP(curblock.block(i,j).data()->domainMap()));
                            break;
                        }
                    }
                }
            }
        }
    }

    for (auto r : ranges)
    {
        if (!r)
            throw new Exception("Error in BlockMatrix::finalize: empty blocks!");
    }

    for (auto d : domains)
    {
        if (!d)
            throw new Exception("Error in BlockMatrix::finalize: empty blocks!");
    }

    for (unsigned int i = 0; i < M_nRows; i++)
    {
        for (unsigned int j = 0; j < M_nCols; j++)
        {
            if (!block(i,j).data())
            {
                block(i,j).data().reset(new MATRIXEPETRA(*ranges[i]));
                block(i,j).data()->zero();
                block(i,j).data()->globalAssemble(domains[j], ranges[i]);
            }
        }
    }

    M_isNull = true;
    for (unsigned int i = 0; i < M_nRows; i++)
    {
        for (unsigned int j = 0; j < M_nCols; j++)
        {
            M_isNull = M_isNull && block(i,j).isNull();
        }
    }

    M_isFinalized = true;
}

template <>
void
BlockMatrix<BlockMatrix<MatrixEp>>::
finalize(BlockMatrix<BlockMatrix<MatrixEp>>* father,
         unsigned int* myRow, unsigned int* myCol)
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

    for (auto r : rows)
    {
        if (r == 0)
            throw new Exception("Error in BlockMatrix::finalize: empty blocks!");
    }

    for (auto c : cols)
    {
        if (c == 0)
            throw new Exception("Error in BlockMatrix::finalize: empty blocks!");
    }

    for (unsigned int i = 0; i < M_nRows; i++)
    {
        for (unsigned int j = 0; j < M_nCols; j++)
        {
            if (block(i,j).nRows() == 0)
            {
                block(i,j).resize(rows[i],rows[j]);
            }
            block(i,j).finalize(this, &i, &j);
        }
    }

    M_isNull = true;
    for (unsigned int i = 0; i < M_nRows; i++)
    {
        for (unsigned int j = 0; j < M_nCols; j++)
        {
            M_isNull = M_isNull && block(i,j).isNull();
        }
    }

    M_isFinalized = true;
}

template <>
void
BlockMatrix<BlockMatrix<MatrixEp>>::
printPattern() const
{
    for (unsigned int i = 0; i < M_nRows; i++)
    {
        for (unsigned int j = 0; j < M_nCols; j++)
        {
            auto curblock = block(i,j);
            if (!curblock.nRows() == 0)
                std::cout << "(" << curblock.nRows() << "," << curblock.nCols() << ")";
            else
                std::cout << "o";
            std::cout << "\t\t";
        }
        std::cout << "\n";
    }
}

template <>
void
BlockMatrix<MatrixEp>::
printPattern() const
{
    for (unsigned int i = 0; i < M_nRows; i++)
    {
        for (unsigned int j = 0; j < M_nCols; j++)
        {
            auto curblock = block(i,j);
            if (curblock.data())
                std::cout << "(" << curblock.data()->rangeMap().mapSize() <<
                             "," << curblock.data()->domainMap().mapSize() << ")";
            else
                std::cout << "o";
            std::cout << "\t\t";
        }
        std::cout << "\n";
    }
}

}
