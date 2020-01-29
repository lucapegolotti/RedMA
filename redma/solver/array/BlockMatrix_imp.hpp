// implementation of BlockMatrix.hpp

namespace RedMA
{

template <class InMatrixType>
BlockMatrix<InMatrixType>::
BlockMatrix() :
  M_nRows(0),
  M_nCols(0),
  M_isFinalized(false),
  M_isNull(false)
{
}

template <class InMatrixType>
BlockMatrix<InMatrixType>::
BlockMatrix(const BlockMatrix& other)
{
    softCopy(other);
}

template <class InMatrixType>
BlockMatrix<InMatrixType>::
BlockMatrix(const unsigned int& nRows, const unsigned int& nCols) :
  M_nRows(nRows),
  M_nCols(nCols)
{
    resize(nRows, nCols);
}

template <class InMatrixType>
BlockMatrix<InMatrixType>
BlockMatrix<InMatrixType>::
getSubmatrix(const unsigned int& ibegin, const unsigned int& iend,
             const unsigned int& jbegin, const unsigned int& jend) const
{
    if (!isFinalized())
        throw new Exception("Matrix must be finalized to get submatrix!");

    BlockMatrix<InMatrixType> retMatrix;

    unsigned int nrows = iend-ibegin+1;
    unsigned int ncols = jend-jbegin+1;
    retMatrix.resize(iend-ibegin+1, jend-jbegin+1);

    for (unsigned int i = ibegin; i <= iend; i++)
    {
        for (unsigned int j = jbegin; j <= jend; j++)
        {
            retMatrix.block(i-ibegin,j-jbegin).softCopy(block(i,j));
        }
    }

    retMatrix.M_isFinalized = true;

    return retMatrix;
}

template <class InMatrixType>
void
BlockMatrix<InMatrixType>::
resize(const unsigned int& nRows, const unsigned int& nCols)
{
    M_nRows = nRows;
    M_nCols = nCols;

    M_matrixGrid.resize(nRows,nCols);
}

template <class InMatrixType>
BlockMatrix<InMatrixType>
BlockMatrix<InMatrixType>::
operator+(const BlockMatrix<InMatrixType>& other) const
{
    if (M_nRows != other.M_nRows || M_nCols != other.M_nCols)
    {
        throw new Exception("Dimension of matrices being added is not consistent!");
    }

    BlockMatrix<InMatrixType> ret(M_nRows,M_nCols);
    ret.hardCopy(*this);

    ret += other;

    return ret;
}

template <class InMatrixType>
void
BlockMatrix<InMatrixType>::
hardCopy(const BlockMatrix<InMatrixType>& other)
{
    resize(other.M_nRows,other.M_nCols);

    for (unsigned int i = 0; i < M_nRows; i++)
    {
        for (unsigned int j = 0; j < M_nCols; j++)
        {
            block(i,j).hardCopy(other.block(i,j));
        }
    }
    M_isFinalized = other.M_isFinalized;
    M_isNull = other.M_isNull;
}

template <class InMatrixType>
void
BlockMatrix<InMatrixType>::
softCopy(const BlockMatrix<InMatrixType>& other)
{
    resize(other.M_nRows,other.M_nCols);

    for (unsigned int i = 0; i < M_nRows; i++)
    {
        for (unsigned int j = 0; j < M_nCols; j++)
        {
            block(i,j).softCopy(other.block(i,j));
        }
    }
    M_isFinalized = other.M_isFinalized;
    M_isNull = other.M_isNull;
}

template <class InMatrixType>
bool
BlockMatrix<InMatrixType>::
isNull() const
{
    if (!M_isFinalized)
        throw new Exception("Matrix must be finalized in order to check if null");
    return M_isNull;
}

template <class InMatrixType>
InMatrixType&
BlockMatrix<InMatrixType>::
block(const unsigned int& iblock, const unsigned int& jblock)
{
    return M_matrixGrid(iblock,jblock);
}

template <class InMatrixType>
InMatrixType
BlockMatrix<InMatrixType>::
block(const unsigned int& iblock, const unsigned int& jblock) const
{
    return M_matrixGrid(iblock,jblock);
}

template <class InMatrixType>
void
BlockMatrix<InMatrixType>::
sumMatrix(const BlockMatrix<InMatrixType>& other)
{
    for (unsigned int i = 0; i < M_nRows; i++)
    {
        for (unsigned int j = 0; j < M_nCols; j++)
        {
            block(i,j) += other.block(i,j);
        }
    }
}

template <class InMatrixType>
BlockMatrix<InMatrixType>&
BlockMatrix<InMatrixType>::
operator+=(const BlockMatrix<InMatrixType>& other)
{
    if (M_nRows == 0 && M_nCols == 0)
    {
        hardCopy(other);
        return *this;
    }

    if (other.nRows() == 0 && other.nCols() == 0)
    {
        return *this;
    }

    if (M_nRows != other.M_nRows || M_nCols != other.M_nCols)
    {
        throw new Exception("Dimension of matrices being added is not consistent!");
    }

    sumMatrix(other);
    return *this;
}

template <class InMatrixType>
void
BlockMatrix<InMatrixType>::
multiplyCoeff(const double& coeff)
{
    for (unsigned int i = 0; i < M_nRows; i++)
    {
        for (unsigned int j = 0; j < M_nCols; j++)
        {
            block(i,j) *= coeff;
        }
    }
}

template <class InMatrixType>
BlockMatrix<InMatrixType>&
BlockMatrix<InMatrixType>::
operator*=(const double& coeff)
{
    multiplyCoeff(coeff);

    return *this;
}

template <class InMatrixType>
template <class InVectorType>
BlockVector<InVectorType>
BlockMatrix<InMatrixType>::
operator*(const BlockVector<InVectorType>& vector) const
{
    if (M_nCols == 0)
        return BlockVector<InVectorType>();

    if (M_nCols != vector.nRows())
        throw new Exception("Dimensions of matrix-vector multiplication are not consistent!");

    BlockVector<InVectorType> ret(M_nRows);

    for (unsigned int i = 0; i < M_nRows; i++)
    {
        for (unsigned int j = 0; j < M_nCols; j++)
        {
            ret.block(i) += block(i,j) * vector.block(j);
        }
    }

    return ret;
}

template <class InMatrixType>
BlockMatrix<InMatrixType>
BlockMatrix<InMatrixType>::
operator*(const BlockMatrix<InMatrixType>& matrix) const
{
    BlockMatrix<InMatrixType> retMatrix;
    if (isNull() || matrix.isNull())
    {
        return retMatrix;
    }

    retMatrix.resize(M_nRows, matrix.nCols());

    for (unsigned int i = 0; i < M_nRows; i++)
    {
        for (unsigned int j = 0; j < matrix.nCols(); j++)
        {
            InMatrixType& curMatrix = retMatrix.block(i,j);
            for (unsigned int k = 0; k < M_nCols; k++)
            {
                curMatrix += block(i,k) * matrix.block(k,j);
            }
        }
    }
    retMatrix.finalize();

    return retMatrix;
}


template <class InMatrixType>
BlockMatrix<InMatrixType>
BlockMatrix<InMatrixType>::
operator*(const double& coeff) const
{
    BlockMatrix<InMatrixType> ret(M_nRows,M_nCols);
    ret.hardCopy(*this);

    ret *= coeff;

    return ret;
}

template <class InMatrixType>
BlockMatrix<InMatrixType>&
BlockMatrix<InMatrixType>::
operator=(const BlockMatrix<InMatrixType>& other)
{
    hardCopy(other);
    return *this;
}

}
