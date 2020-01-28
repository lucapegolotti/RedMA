// implementation of BlockMatrix.hpp

namespace RedMA
{

template <class InMatrixType>
BlockMatrix<InMatrixType>::
BlockMatrix() :
  M_nRows(0),
  M_nCols(0),
  M_isFinalized(true)
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


    for (unsigned int i = 0; i < M_nRows; i++)
    {
        for (unsigned int j = 0; j < M_nCols; j++)
        {
            block(i,j) += other.block(i,j);
        }
    }
}

template <class InMatrixType>
void
BlockMatrix<InMatrixType>::
printPattern() const
{
    for (unsigned int i = 0; i < M_nRows; i++)
    {
        for (unsigned int j = 0; j < M_nCols; j++)
        {
            auto curblock = block(i,j);
            if (!curblock.isNull())
                std::cout << "x";
            else
                std::cout << "o";
            std::cout << "\t";
        }
        std::cout << "\n";
    }
}

template <class InMatrixType>
BlockMatrix<InMatrixType>&
BlockMatrix<InMatrixType>::
operator*=(const double& coeff)
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

template <class InMatrixType>
void
BlockMatrix<InMatrixType>::
finalize()
{
    M_isFinalized = true;
}

}
