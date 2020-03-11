#include "BlockVector.hpp"

namespace RedMA
{

template <class InVectorType>
BlockVector<InVectorType>::
BlockVector() :
  M_nRows(0)
{
}

template <class InVectorType>
BlockVector<InVectorType>::
BlockVector(BlockVector& other)
{
    softCopy(other);
}

template <class InVectorType>
BlockVector<InVectorType>::
BlockVector(const BlockVector& other)
{
    softCopy(other);
}

template <class InVectorType>
BlockVector<InVectorType>::
BlockVector(const unsigned int& nRows) :
  M_nRows(nRows)
{
    resize(nRows);
}

template <class InVectorType>
void
BlockVector<InVectorType>::
resize(const unsigned int& nRows)
{
    M_nRows = nRows;
    M_vectorGrid.resize(nRows,1);
}

template <class InVectorType>
BlockVector<InVectorType>
BlockVector<InVectorType>::
operator*(const double& coeff) const
{
    BlockVector<InVectorType> ret;
    ret.hardCopy(*this);

    ret *= coeff;

    return ret;
}

template <class InVectorType>
double
BlockVector<InVectorType>::
norm2() const
{
    double ret = 0;

    for (unsigned int i = 0; i < M_nRows; i++)
    {
        double curNorm = block(i).norm2();
        ret += curNorm * curNorm;
    }

    return sqrt(ret);
}

template <class InVectorType>
BlockVector<InVectorType>&
BlockVector<InVectorType>::
operator*=(const double& coeff)
{
    for (unsigned int i = 0; i < M_nRows; i++)
        block(i) *= coeff;

    return *this;
}

template <class InVectorType>
BlockVector<InVectorType>
BlockVector<InVectorType>::
operator+(const BlockVector<InVectorType>& other) const
{
    BlockVector<InVectorType> ret;

    if (M_nRows == 0)
    {
        ret.hardCopy(other);
        return ret;
    }

    if (other.nRows() == 0)
    {
        ret.hardCopy(*this);
        return ret;
    }

    if (M_nRows != other.M_nRows)
    {
        throw new Exception("Dimension of vectors being added is not consistent!");
    }

    ret.hardCopy(*this);

    ret += other;
    return ret;
}

template <class InVectorType>
BlockVector<InVectorType>
BlockVector<InVectorType>::
operator-(const BlockVector<InVectorType>& other) const
{
    BlockVector<InVectorType> ret;

    if (M_nRows == 0)
    {
        ret.hardCopy(other);
        return ret;
    }

    if (other.nRows() == 0)
    {
        ret.hardCopy(*this);
        return ret;
    }

    if (M_nRows != other.M_nRows)
    {
        throw new Exception("Dimension of vectors being subtracted is not consistent!");
    }

    ret.hardCopy(*this);

    ret -= other;
    return ret;
}

template <class InVectorType>
BlockVector<InVectorType>&
BlockVector<InVectorType>::
operator+=(const BlockVector<InVectorType>& other)
{
    if (M_nRows == 0)
    {
        hardCopy(other);
        return *this;
    }

    if (other.nRows() == 0)
        return *this;

    if (M_nRows != other.M_nRows)
    {
        throw new Exception("Dimension of vectors being added is not consistent!");
    }

    for (unsigned int i = 0; i < M_nRows; i++)
        block(i) += other.block(i);

    return *this;
}

template <class InVectorType>
BlockVector<InVectorType>&
BlockVector<InVectorType>::
operator-=(const BlockVector<InVectorType>& other)
{
    if (M_nRows == 0)
    {
        hardCopy(other);
        (*this) *= (-1);
        return *this;
    }

    if (other.nRows() == 0)
        return *this;

    if (M_nRows != other.M_nRows)
    {
        throw new Exception("Dimension of vectors being subtracted is not consistent!");
    }

    for (unsigned int i = 0; i < M_nRows; i++)
        block(i) -= other.block(i);

    return *this;
}

template <class InVectorType>
void
BlockVector<InVectorType>::
zero()
{
    for (unsigned int i = 0; i < M_nRows; i++)
        block(i) *= 0;
}

template <class InVectorType>
void
BlockVector<InVectorType>::
hardCopy(const BlockVector<InVectorType>& other)
{
    resize(other.M_nRows);

    for (unsigned int i = 0; i < M_nRows; i++)
        block(i).hardCopy(other.block(i));
}

template <class InVectorType>
void
BlockVector<InVectorType>::
softCopy(const BlockVector<InVectorType>& other)
{
    resize(other.M_nRows);
    for (unsigned int i = 0; i < M_nRows; i++)
    {
        block(i).softCopy(other.block(i));
    }
}

template <class InVectorType>
InVectorType&
BlockVector<InVectorType>::
block(const unsigned int& iblock)
{
    return M_vectorGrid(iblock,0);
}

template <class InVectorType>
InVectorType
BlockVector<InVectorType>::
block(const unsigned int& iblock) const
{
    return M_vectorGrid(iblock,0);
}

template <class InVectorType>
BlockVector<InVectorType>
BlockVector<InVectorType>::
getSubvector(const unsigned int& ibegin, const unsigned int& iend) const
{
    BlockVector<InVectorType> retVector;

    unsigned int nrows = iend-ibegin+1;
    retVector.resize(iend-ibegin+1);

    for (unsigned int i = ibegin; i <= iend; i++)
        retVector.block(i-ibegin).softCopy(block(i));

    return retVector;
}

}
