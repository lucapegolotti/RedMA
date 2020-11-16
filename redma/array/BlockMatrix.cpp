#include "BlockMatrix.hpp"

namespace RedMA
{

BlockMatrix::
BlockMatrix() :
  aMatrix(BLOCK),
  M_isOpen(true)
{

}

BlockMatrix::
BlockMatrix(const BlockMatrix& other) :
  aMatrix(BLOCK),
  M_isOpen(true)
{
    resize(other.nRows(), other.nCols());

    for (unsigned int i = 0; i < other.nRows(); i++)
    {
        for (unsigned int j = 0; j < other.nCols(); j++)
        {
            setBlock(i,j,other.block(i,j));
        }
    }
}

BlockMatrix::
BlockMatrix(const unsigned int& nRows, const unsigned int& nCols) :
  aMatrix(BLOCK),
  M_isOpen(true)
{
    resize(nRows, nCols);
}

void
BlockMatrix::
add(std::shared_ptr<aMatrix> other)
{
    std::shared_ptr<BlockMatrix> otherMatrix = std::static_pointer_cast<BlockMatrix>(other);

    if (other->isZero())
        return;

    if (nRows() != other->nRows() || nCols() != other->nCols())
        throw new Exception("BlockMatrix: inconsistent dimensions in add!");

    for (unsigned int i = 0; i < nRows(); i++)
    {
        for (unsigned int j = 0; j < nCols(); j++)
        {
            if (!block(i,j)->isZero())
                block(i,j)->add(otherMatrix->block(i,j));
            else
            {
                if (other->block(i,j))
                    setBlock(i,j,std::shared_ptr<aMatrix>(other->block(i,j)->clone()));
            }
        }
    }
    updateNormInf();
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
    updateNormInf();
}

std::shared_ptr<aMatrix>
BlockMatrix::
multiplyByMatrix(std::shared_ptr<aMatrix> other)
{
    std::shared_ptr<BlockMatrix> otherMatrix = std::static_pointer_cast<BlockMatrix>(other);
    std::shared_ptr<BlockMatrix> retMatrix(new BlockMatrix(nRows(), other->nCols()));

    if (isZero() || other->isZero())
    {
        retMatrix.reset(new BlockMatrix(0, 0));
        return retMatrix;
    }

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
    retMatrix->close();
    retMatrix->updateNormInf();
    return retMatrix;
}

std::shared_ptr<aMatrix>
BlockMatrix::
transpose() const
{
    // checkClosed();

    std::shared_ptr<BlockMatrix> retMatrix(new BlockMatrix(nCols(), nRows()));

    for (unsigned int i = 0; i < nRows(); i++)
    {
        for (unsigned int j = 0; j < nCols(); j++)
        {
            retMatrix->setBlock(j,i,block(i,j)->transpose());
        }
    }

    retMatrix->close();

    return retMatrix;
}

std::shared_ptr<aVector>
BlockMatrix::
multiplyByVector(std::shared_ptr<aVector> vector)
{
    checkType(vector, BLOCK);

    BlockVector* otherVector = dynamic_cast<BlockVector*>(vector.get());
    // checkClosed();
    // otherVector->checkClosed();

    std::shared_ptr<BlockVector> retVector(new BlockVector(nRows()));

    for (unsigned int i = 0; i < nRows(); i++)
    {
        for (unsigned int j = 0; j < nCols(); j++)
        {
            // check if we need to the conversion distributed->dense
            std::shared_ptr<aVector> tempRes;
            if (!block(i,j)->isZero())
            {
                if (block(i,j)->type() == DENSE && otherVector->block(j)->type() == DISTRIBUTED)
                {
                    // printlog(YELLOW, "[BlockMatrix::multiplyByVector] explicit dense->distributed conversion\n", true);
                    std::shared_ptr<DistributedVector> asDistributed = std::static_pointer_cast<DistributedVector>(otherVector->block(j));
                    std::shared_ptr<DenseVector> asDense = asDistributed->toDenseVectorPtr();
                    tempRes = DistributedVector::convertDenseVector(std::static_pointer_cast<DenseVector>(block(i,j)->multiplyByVector(asDense)),asDistributed->commPtr());
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
    // retVector->close();

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
                    if (!std::static_pointer_cast<BlockMatrix>(block(i,j))->globalTypeIs(type))
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
                    M_comm = std::static_pointer_cast<BlockMatrix>(block(i,j))->commPtr();
                else if (!block(i,j)->isZero() && block(i,j)->type() == SPARSE)
                    M_comm = std::static_pointer_cast<SparseMatrix>(block(i,j))->commPtr();
            }
        }
    }
}

std::shared_ptr<BlockMatrix>
BlockMatrix::
convertInnerTo(Datatype type, std::shared_ptr<Epetra_Comm> comm)
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

    std::shared_ptr<BlockMatrix> retMat(new BlockMatrix(nRows(),nCols()));

    for (unsigned int i = 0; i < nRows(); i++)
    {
        for (unsigned int j = 0; j < nCols(); j++)
        {
            retMat->setBlock(i,j,block(i,j));
            if (block(i,j))
            {
                if (!block(i,j)->isZero() && block(i,j)->type() == BLOCK)
                {
                    retMat->setBlock(i,j,std::static_pointer_cast<BlockMatrix>(block(i,j))->convertInnerTo(type,comm));
                }
                else if (!block(i,j)->isZero() && block(i,j)->type() != type)
                {
                    if (block(i,j)->type() != DENSE)
                        throw new Exception("All inner matrices which are not sparse must be dense!");
                    retMat->setBlock(i,j,SparseMatrix::convertDenseMatrix(std::static_pointer_cast<DenseMatrix>(block(i,j)),comm));
                }
            }
        }
    }
    retMat->close();

    return retMat;
}

void
BlockMatrix::
softCopy(std::shared_ptr<aMatrix> other)
{
    if (other)
    {
        checkType(other, BLOCK);
        auto otherMatrix = dynamic_cast<BlockMatrix*>(other.get());

        resize(other->nRows(), other->nCols());

        for (unsigned int i = 0; i < nRows(); i++)
        {
            for (unsigned int j = 0; j < nCols(); j++)
            {
                M_matrixGrid(i,j)->softCopy(otherMatrix->block(i,j));
            }
        }

        M_isOpen = otherMatrix->isOpen();
        M_normInf = otherMatrix->normInf();
    }
}

void
BlockMatrix::
hardCopy(std::shared_ptr<aMatrix> other)
{
    if (other)
    {
        checkType(other, BLOCK);
        auto otherMatrix = dynamic_cast<BlockMatrix*>(other.get());

        resize(other->nRows(), other->nCols());

        for (unsigned int i = 0; i < nRows(); i++)
        {
            for (unsigned int j = 0; j < nCols(); j++)
            {
                if (M_matrixGrid(i,j)->isZero())
                    M_matrixGrid(i,j).reset(otherMatrix->block(i,j)->clone());
                else
                    M_matrixGrid(i,j)->hardCopy(otherMatrix->block(i,j));
            }
        }

        M_isOpen = otherMatrix->isOpen();
        M_normInf = otherMatrix->normInf();
    }
}

bool
BlockMatrix::
isZero()
{
    for (unsigned int i = 0; i < nRows(); i++)
    {
        for (unsigned int j = 0; j < nCols(); j++)
        {
            if (block(i,j) && !block(i,j)->isZero())
                return false;
        }
    }
    return M_normInf < ZEROTHRESHOLD;
}

void
BlockMatrix::
dump(std::string filename) const
{
    throw new Exception("Dump not implemented for BlockMatrix");
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
            retMatrix->setBlock(i,j,std::shared_ptr<aMatrix>(M_matrixGrid(i,j)->clone()));
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

std::shared_ptr<aMatrix>
BlockMatrix::
block(const unsigned int& iblock, const unsigned int& jblock) const
{
    return M_matrixGrid(iblock,jblock);
}

std::shared_ptr<BlockMatrix>
BlockMatrix::
getSubmatrix(const unsigned int& ibegin, const unsigned int& iend,
             const unsigned int& jbegin, const unsigned int& jend) const
{
    std::shared_ptr<BlockMatrix> retMatrix(new BlockMatrix());

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

    retMatrix->close();

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
         std::shared_ptr<aMatrix> matrix)
{
    if (matrix)
    {
        checkOpen();
        M_matrixGrid(iblock,jblock) = matrix;
        updateNormInf();
    }
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
                return std::static_pointer_cast<BlockMatrix>(block(i,j))->level() + 1;
        }
    }
    return 1;
}


void
BlockMatrix::
updateNormInf()
{
    for (unsigned int i = 0; i < nRows(); i++)
    {
        for (unsigned int j = 0; j < nCols(); j++)
        {
            if (block(i,j))
                M_normInf = M_normInf < block(i,j)->normInf() ? block(i,j)->normInf() : M_normInf;
        }
    }
}

void
BlockMatrix::
close()
{
    M_isOpen = false;

    for (unsigned int i = 0; i < nRows(); i++)
    {
        for (unsigned int j = 0; j < nCols(); j++)
        {
            if (block(i,j) == nullptr)
                block(i,j).reset(new DenseMatrix());
        }
    }
}

}
