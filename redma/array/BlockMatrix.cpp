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
    // other.checkClosed();

    resize(other.nRows(), other.nCols());

    for (unsigned int i = 0; i < other.nRows(); i++)
    {
        for (unsigned int j = 0; j < other.nCols(); j++)
        {
            setBlock(i,j,other.block(i,j));
        }
    }
    // close();
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

    // checkClosed();
    // otherMatrix->checkClosed();

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
    // checkClosed();

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
    // checkClosed();
    // otherMatrix->checkClosed();
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
            // std::shared_ptr<aMatrix>& curMatrix = retMatrix->block(i,j);
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
            if (block(i,j)->type() == DENSE && otherVector->block(j)->type() == DISTRIBUTED)
            {
                printlog(YELLOW, "[BlockMatrix::multiplyByVector] explicit dense->distributed conversion\n", true);
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
    // checkClosed();

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

// BlockMatrix::
// BlockMatrix() :
//   M_nRows(0),
//   M_nCols(0),
//   M_isFinalized(false),
//   M_isNull(false)
// {
// }
//
// BlockMatrix::
// BlockMatrix(const BlockMatrix& other)
// {
//     softCopy(other);
// }
//
// BlockMatrix::
// BlockMatrix(const unsigned int& nRows, const unsigned int& nCols) :
//   M_nRows(nRows),
//   M_nCols(nCols)
// {
//     resize(nRows, nCols);
// }
//
// BlockMatrix
// BlockMatrix::
// getSubmatrix(const unsigned int& ibegin, const unsigned int& iend,
//              const unsigned int& jbegin, const unsigned int& jend) const
// {
//     if (!isFinalized())
//         throw new Exception("Matrix must be finalized to get submatrix!");
//
//     BlockMatrix retMatrix;
//
//     unsigned int nrows = iend-ibegin+1;
//     unsigned int ncols = jend-jbegin+1;
//     retMatrix.resize(iend-ibegin+1, jend-jbegin+1);
//
//     for (unsigned int i = ibegin; i <= iend; i++)
//     {
//         for (unsigned int j = jbegin; j <= jend; j++)
//         {
//             retMatrix.block(i-ibegin,j-jbegin)->softCopy(block(i,j));
//         }
//     }
//
//     retMatrix.finalize();
//
//     return retMatrix;
// }
//
// void
// BlockMatrix::
// resize(const unsigned int& nRows, const unsigned int& nCols)
// {
//     M_nRows = nRows;
//     M_nCols = nCols;
//
//     M_matrixGrid.resize(nRows,nCols);
// }
//
// BlockMatrix
// BlockMatrix::
// operator+(const BlockMatrix& other) const
// {
//     if (M_nRows != other.M_nRows || M_nCols != other.M_nCols)
//         throw new Exception("Dimension of matrices being added is not consistent!");
//
//     BlockMatrix ret(M_nRows,M_nCols);
//     ret.hardCopy(*this);
//
//     ret += other;
//
//     return ret;
// }
//
// void
// BlockMatrix::
// hardCopy(const BlockMatrix& other)
// {
//     resize(other.M_nRows,other.M_nCols);
//
//     for (unsigned int i = 0; i < M_nRows; i++)
//     {
//         for (unsigned int j = 0; j < M_nCols; j++)
//         {
//             block(i,j)->hardCopy(other.block(i,j));
//         }
//     }
//     M_isFinalized = other.M_isFinalized;
//     M_isNull = other.M_isNull;
// }
//
// void
// BlockMatrix::
// softCopy(const BlockMatrix& other)
// {
//     resize(other.M_nRows,other.M_nCols);
//
//     for (unsigned int i = 0; i < M_nRows; i++)
//     {
//         for (unsigned int j = 0; j < M_nCols; j++)
//         {
//             block(i,j)->softCopy(other.block(i,j));
//         }
//     }
//     M_isFinalized = other.M_isFinalized;
//     M_isNull = other.M_isNull;
// }
//
// bool
// BlockMatrix::
// isNull() const
// {
//     if (!M_isFinalized)
//         throw new Exception("Matrix must be finalized in order to check if null");
//     return M_isNull;
// }
//
// std::shared_ptr<aMatrix>&
// BlockMatrix::
// block(const unsigned int& iblock, const unsigned int& jblock)
// {
//     return M_matrixGrid(iblock,jblock);
// }
//
// std::shared_ptr<aMatrix>
// BlockMatrix::
// block(const unsigned int& iblock, const unsigned int& jblock) const
// {
//     return M_matrixGrid(iblock,jblock);
// }
//
// void
// BlockMatrix::
// sumMatrix(const BlockMatrix& other)
// {
//     for (unsigned int i = 0; i < M_nRows; i++)
//     {
//         for (unsigned int j = 0; j < M_nCols; j++)
//         {
//             block(i,j)->add(other.block(i,j));
//         }
//     }
// }
//
// BlockMatrix&
// BlockMatrix::
// operator+=(const BlockMatrix& other)
// {
//     if (M_nRows == 0 && M_nCols == 0)
//     {
//         hardCopy(other);
//         return *this;
//     }
//
//     if (other.nRows() == 0 && other.nCols() == 0)
//         return *this;
//
//     if (M_nRows != other.M_nRows || M_nCols != other.M_nCols)
//         throw new Exception("Dimension of matrices being added is not consistent!");
//
//     sumMatrix(other);
//     return *this;
// }

// BlockMatrix&
// BlockMatrix::
// operator-=(const BlockMatrix& other)
// {
//     if (M_nRows == 0 && M_nCols == 0)
//     {
//         hardCopy(other);
//         *this *= (-1.0);
//         return *this;
//     }
//
//     if (other.nRows() == 0 && other.nCols() == 0)
//         return *this;
//
//     if (M_nRows != other.M_nRows || M_nCols != other.M_nCols)
//         throw new Exception("Dimension of matrices being subtracted is not consistent!");
//
//     subtractMatrix(other);
//     return *this;
// }

// void
// BlockMatrix::
// multiplyCoeff(const double& coeff)
// {
//     for (unsigned int i = 0; i < M_nRows; i++)
//     {
//         for (unsigned int j = 0; j < M_nCols; j++)
//         {
//             block(i,j)->multiplyByScalar(coeff);
//         }
//     }
// }
//
// BlockMatrix&
// BlockMatrix::
// operator*=(const double& coeff)
// {
//     multiplyCoeff(coeff);
//
//     return *this;
// }
//
// BlockVector
// BlockMatrix::
// operator*(const BlockVector& vector) const
// {
//     if (M_nCols == 0)
//         return BlockVector();
//
//     if (M_nCols != vector.nRows())
//         throw new Exception("Dimensions of matrix-vector multiplication are not consistent!");
//
//     BlockVector ret(M_nRows);
//
//     for (unsigned int i = 0; i < M_nRows; i++)
//     {
//         for (unsigned int j = 0; j < M_nCols; j++)
//         {
//             aMatrix* copy = block(i,j)->clone();
//             std::shared_ptr<aVector> tempVec = copy->multiplyByVector(vector.block(j));
//             ret.block(i)->add(tempVec);
//             // the memory management is left to us
//             delete copy;
//         }
//     }
//
//     return ret;
// }
//
// BlockMatrix
// BlockMatrix::
// operator*(const BlockMatrix& matrix) const
// {
//     BlockMatrix retMatrix;
//
//     if (isNull() || matrix.isNull())
//         return retMatrix;
//
//     retMatrix.resize(M_nRows, matrix.nCols());
//
//     for (unsigned int i = 0; i < M_nRows; i++)
//     {
//         for (unsigned int j = 0; j < matrix.nCols(); j++)
//         {
//             std::shared_ptr<aMatrix> curMatrix = retMatrix.block(i,j);
//             for (unsigned int k = 0; k < M_nCols; k++)
//             {
//                 aMatrix* copy = block(i,k)->clone();
//                 curMatrix->add(copy->multiplyByMatrix(matrix.block(k,j)));
//                 delete copy;
//             }
//         }
//     }
//     retMatrix.finalize();
//
//     return retMatrix;
// }
//
//
// BlockMatrix
// BlockMatrix::
// operator*(const double& coeff) const
// {
//     BlockMatrix ret(M_nRows,M_nCols);
//     ret.hardCopy(*this);
//
//     ret *= coeff;
//
//     return ret;
// }
//
// BlockMatrix&
// BlockMatrix::
// operator=(const BlockMatrix& other)
// {
//     hardCopy(other);
//     return *this;
// }

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
    // M_dimensionsRows.resize(nRows());
    // M_dimensionsCols.resize(nCols());
    // // build block dimensions
    // for (unsigned int i = 0; i < nRows(); i++)
    // {
    //     for (unsigned int j = 0; j < nCols(); j++)
    //     {
    //         if (!block(i,j)->isZero() && block(i,j)->type() == BLOCK)
    //         {
    //             SHP(BlockMatrix) curMatrix = std::static_pointer_cast<BlockMatrix>(block(i,j));
    //             curMatrix->close();
    //             M_dimensionsRows[i] = curMatrix->getDimensionsRows();
    //         }
    //     }
    // }
}

}
