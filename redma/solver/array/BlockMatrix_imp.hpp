// implementation of BlockMatrix.hpp

namespace RedMA
{

template <class InMatrixType>
BlockMatrix<InMatrixType>::
BlockMatrix() :
  M_nRows(0),
  M_nCols(0)
{
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
template <class VectorType>
VectorType
BlockMatrix<InMatrixType>::
operator*(const VectorType& vector)
{

}

}
