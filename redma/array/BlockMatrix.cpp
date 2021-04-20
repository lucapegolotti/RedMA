#include "BlockMatrix.hpp"

namespace RedMA
{

BlockMatrix::
BlockMatrix()
{

}


BlockMatrix::
BlockMatrix(const unsigned int& nRows, const unsigned int& nCols)
{
    resize(nRows, nCols);
}

void
BlockMatrix::
add(shp<aMatrix> other)
{
    if (other->isZero())
        return;

    if (isZero())
    {
        deepCopy(other);
        return;
    }

    if (nRows() != other->nRows() || nCols() != other->nCols())
        throw new Exception("BlockMatrix: inconsistent dimensions in add!");

    shp<BlockMatrix> otherMatrix = convert<BlockMatrix>(other);

    for (unsigned int i = 0; i < nRows(); i++)
    {
        for (unsigned int j = 0; j < nCols(); j++)
        {
            if (!block(i,j)->isZero())
                block(i,j)->add(otherMatrix->block(i,j));
            else
            {
                if (otherMatrix->block(i,j))
                {
                    shp<aDataWrapper> blck(otherMatrix->block(i,j)->clone());
                    setBlock(i,j,spcast<aMatrix>(blck));
                }
            }
        }
    }
}

void
BlockMatrix::
multiplyByScalar(const double& coeff)
{
    if (isZero())
        return;

    for (int i = 0; i < nRows(); i++)
    {
        for (int j = 0; j < nCols(); j++)
        {
            if (block(i,j))
                block(i,j)->multiplyByScalar(coeff);
        }
    }
}

shp<aMatrix>
BlockMatrix::
multiplyByMatrix(shp<aMatrix> other)
{
    shp<BlockMatrix> retMatrix(new BlockMatrix(nRows(), other->nCols()));

    if (isZero() || other->isZero())
    {
        retMatrix.reset(new BlockMatrix(0, 0));
        return retMatrix;
    }

    shp<BlockMatrix> otherMatrix = convert<BlockMatrix>(other);

    for (unsigned int i = 0; i < nRows(); i++)
    {
        for (unsigned int j = 0; j < other->nCols(); j++)
        {
            for (unsigned int k = 0; k < nCols(); k++)
            {
                if (block(i,k))
                {
                    if (retMatrix->block(i,j)->isZero())
                        retMatrix->setBlock(i,j,block(i,k)->multiplyByMatrix(otherMatrix->block(k,j)));
                    else
                        retMatrix->block(i,j)->add(block(i,k)->multiplyByMatrix(otherMatrix->block(k,j)));
                }
            }
        }
    }
    return retMatrix;
}

shp<aMatrix>
BlockMatrix::
transpose() const
{
    shp<BlockMatrix> retMatrix(new BlockMatrix(nCols(), nRows()));

    for (unsigned int i = 0; i < nRows(); i++)
    {
        for (unsigned int j = 0; j < nCols(); j++)
        {
            retMatrix->setBlock(j,i,block(i,j)->transpose());
        }
    }

    return retMatrix;
}

shp<aVector>
BlockMatrix::
multiplyByVector(shp<aVector> vector)
{
    shp<BlockVector> retVector(new BlockVector(nRows()));

    if (vector->isZero())
        return retVector;

    shp<BlockVector> otherVector = convert<BlockVector>(vector);

    for (unsigned int i = 0; i < nRows(); i++)
    {
        for (unsigned int j = 0; j < nCols(); j++)
        {
            // check if we need to make the conversion distributed->dense
            shp<aVector> tempRes;
            if (!block(i,j)->isZero())
            {
                if (block(i,j)->type() == DENSE && otherVector->block(j)->type() == DISTRIBUTED)
                {
                    // printlog(YELLOW, "[BlockMatrix::multiplyByVector] explicit dense->distributed conversion\n", true);
                    shp<DistributedVector> asDistributed = convert<DistributedVector>(otherVector->block(j));
                    shp<DenseVector> asDense = asDistributed->toDenseVectorPtr();
                    // tempRes = DistributedVector::convertDenseVector(spcast<DenseVector>(block(i,j)->multiplyByVector(asDense)),asDistributed->commPtr());
                    tempRes = block(i,j)->multiplyByVector(asDense);
                }
                else
                    tempRes = block(i,j)->multiplyByVector(otherVector->block(j));

                if (retVector->block(i)->isZero())
                    retVector->setBlock(i,tempRes);
                else
                    retVector->block(i)->add(tempRes);
            }
        }
    }

    return retVector;
}

bool
BlockMatrix::
globalTypeIs(Datatype type)
{
    for (unsigned int i = 0; i < nRows(); i++)
    {
        for (unsigned int j = 0; j < nCols(); j++)
        {
            if (block(i,j))
            {
                if (!block(i,j)->isZero() && block(i,j)->type() == BLOCK)
                {
                    if (!spcast<BlockMatrix>(block(i,j))->globalTypeIs(type))
                        return false;
                }
                else if (!block(i,j)->isZero() && block(i,j)->type() != type)
                {
                    return false;
                }
            }
        }
    }
    return true;
}

void
BlockMatrix::
findComm()
{
    M_comm = nullptr;

    for (unsigned int i = 0; i < nRows(); i++)
    {
        for (unsigned int j = 0; j < nCols(); j++)
        {
            if (block(i,j) && M_comm == nullptr)
            {
                if (!block(i,j)->isZero() && block(i,j)->type() == BLOCK)
                    M_comm = spcast<BlockMatrix>(block(i,j))->commPtr();
                else if (!block(i,j)->isZero() && block(i,j)->type() == SPARSE)
                    M_comm = spcast<SparseMatrix>(block(i,j))->commPtr();
            }
        }
    }
}

shp<BlockMatrix>
BlockMatrix::
convertInnerTo(Datatype type, shp<Epetra_Comm> comm)
{
    if (type != SPARSE)
        throw new Exception("convertInnerTo implemented only for sparse output");

    if (comm == nullptr)
    {
        findComm();
        comm = M_comm;
    }

    if (comm == nullptr)
        throw new Exception("No communicator in block matrix!");

    shp<BlockMatrix> retMat(new BlockMatrix(nRows(),nCols()));

    for (unsigned int i = 0; i < nRows(); i++)
    {
        for (unsigned int j = 0; j < nCols(); j++)
        {
            retMat->setBlock(i,j,block(i,j));
            if (block(i,j))
            {
                if (!block(i,j)->isZero() && block(i,j)->type() == BLOCK)
                {
                    retMat->setBlock(i,j,spcast<BlockMatrix>(block(i,j))->convertInnerTo(type,comm));
                }
                else if (!block(i,j)->isZero() && block(i,j)->type() != type)
                {
                    if (block(i,j)->type() != DENSE)
                        throw new Exception("All inner matrices which are not sparse must be dense!");
                    retMat->setBlock(i,j,SparseMatrix::convertDenseMatrix(spcast<DenseMatrix>(block(i,j)),comm));
                }
            }
        }
    }

    return retMat;
}

void
BlockMatrix::
shallowCopy(shp<aDataWrapper> other)
{
    if (other && !other->isZero())
    {
        auto otherMatrix = convert<BlockMatrix>(other);

        resize(otherMatrix->nRows(), otherMatrix->nCols());

        for (unsigned int i = 0; i < nRows(); i++)
        {
            for (unsigned int j = 0; j < nCols(); j++)
            {
                M_matrixGrid(i,j)->shallowCopy(otherMatrix->block(i,j));
            }
        }
    }
}

void
BlockMatrix::
deepCopy(shp<aDataWrapper> other)
{
    if (other && !other->isZero())
    {
        auto otherMatrix = convert<BlockMatrix>(other);

        resize(otherMatrix->nRows(), otherMatrix->nCols());

        for (unsigned int i = 0; i < nRows(); i++)
        {
            for (unsigned int j = 0; j < nCols(); j++)
            {
                if (M_matrixGrid(i,j)->isZero())
                    M_matrixGrid(i,j).reset(static_cast<aMatrix*>(otherMatrix->block(i,j)->clone()));
                else
                    M_matrixGrid(i,j)->deepCopy(otherMatrix->block(i,j));
            }
        }
    }
}

bool
BlockMatrix::
isZero() const
{
    for (unsigned int i = 0; i < nRows(); i++)
    {
        for (unsigned int j = 0; j < nCols(); j++)
        {
            if (block(i,j) && !block(i,j)->isZero())
                return false;
        }
    }
    return true;
}

void
BlockMatrix::
dump(std::string filename) const
{
    for (unsigned int i = 0; i < nRows(); i++)
    {
        for (unsigned int j = 0; j < nCols(); j++)
        {
            if (!M_matrixGrid(i,j)->isZero())
                M_matrixGrid(i,j)->dump(filename + "_block_" + std::to_string(i) + "_" + std::to_string(j));
        }
    }
}

BlockMatrix*
BlockMatrix::
clone() const
{
    BlockMatrix* retMatrix = new BlockMatrix(nRows(),nCols());
    for (unsigned int i = 0; i < nRows(); i++)
    {
        for (unsigned int j = 0; j < nCols(); j++)
        {
            retMatrix->setBlock(i,j,shp<aMatrix>(static_cast<aMatrix*>(M_matrixGrid(i,j)->clone())));
        }
    }
    return retMatrix;
}

void
BlockMatrix::
resize(const unsigned int& nRows, const unsigned int& nCols)
{
    M_nRows = nRows;
    M_nCols = nCols;
    M_matrixGrid.resize(nRows,nCols);
    for (unsigned int i = 0; i < nRows; i++)
    {
        for (unsigned int j = 0; j < nCols; j++)
        {
            // we use dense matrices as placeholders for zero matrices
            M_matrixGrid(i,j).reset(new DenseMatrix());
        }
    }
}

shp<aMatrix>
BlockMatrix::
block(const unsigned int& iblock, const unsigned int& jblock) const
{
    if (iblock >= nRows())
        throw new Exception("[BlockMatrix::block] iblock >= nRows()");

    if (jblock >= nCols())
        throw new Exception("[BlockMatrix::block] jblock >= nCols()");

    return M_matrixGrid(iblock,jblock);
}

shp<BlockMatrix>
BlockMatrix::
getSubmatrix(const unsigned int& ibegin, const unsigned int& iend,
             const unsigned int& jbegin, const unsigned int& jend) const
{
    if (ibegin > nRows() || iend > nRows() || jbegin > nCols() || jend > nCols())
        throw new Exception("[BlockMatrix::getSubmatrix] wrong bounds!");

    shp<BlockMatrix> retMatrix(new BlockMatrix());

    unsigned int nrows = iend-ibegin+1;
    unsigned int ncols = jend-jbegin+1;
    retMatrix->resize(iend-ibegin+1, jend-jbegin+1);

    for (unsigned int i = ibegin; i <= iend; i++)
    {
        for (unsigned int j = jbegin; j <= jend; j++)
        {
            retMatrix->setBlock(i-ibegin,j-jbegin,block(i,j));
        }
    }


    return retMatrix;
}

void
BlockMatrix::
printPattern() const
{
    for (unsigned int i = 0; i < nRows(); i++)
    {
        for (unsigned int j = 0; j < nCols(); j++)
        {
            if (block(i,j))
            {
                unsigned int nRows = block(i,j)->nRows();
                unsigned int nCols = block(i,j)->nCols();

                std::cout << "(" << nRows << "," << nCols << ")" << "\t";
            }
            else
            {
                std::cout << "(null)" << "\t";
            }
        }
        std::cout << "\n";
    }
}

void
BlockMatrix::
setBlock(const unsigned int& iblock, const unsigned int& jblock,
         shp<aMatrix> matrix)
{
    if (iblock >= nRows())
        throw new Exception("[BlockMatrix::setBlock] iblock >= nRows()");

    if (jblock >= nCols())
        throw new Exception("[BlockMatrix::setBlock] jblock >= nCols()");

    if (matrix)
        M_matrixGrid(iblock,jblock) = matrix;
}

unsigned int
BlockMatrix::
level()
{
    for (unsigned int i = 0; i < nRows(); i++)
    {
        for (unsigned int j = 0; j < nCols(); j++)
        {
            if (block(i,j)->type() == BLOCK)
                return convert<BlockMatrix>(block(i,j))->level() + 1;
        }
    }
    return 1;
}

}
