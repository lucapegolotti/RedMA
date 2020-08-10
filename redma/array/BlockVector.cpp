#include "BlockVector.hpp"

namespace RedMA
{

BlockVector::
BlockVector() :
  aVector(BLOCK),
  M_isOpen(true)
{
}

BlockVector::
BlockVector(const BlockVector& other) :
  aVector(BLOCK),
  M_isOpen(true)
{
    other.checkClosed();

    resize(other.nRows());

    for (unsigned int i = 0; i < other.nRows(); i++)
    {
        setBlock(i, other.block(i));
    }
    close();
}


BlockVector::
BlockVector(const unsigned int& nRows) :
  aVector(BLOCK),
  M_isOpen(true)
{
    resize(nRows);
}

void
BlockVector::
resize(const unsigned int& nRows)
{
    M_nRows = nRows;
    M_vectorGrid.resize(nRows,1);
}

void
BlockVector::
add(std::shared_ptr<aVector> other)
{
    BlockVector* otherVector = dynamic_cast<BlockVector*>(other.get());

    checkClosed();
    otherVector->checkClosed();

    if (nRows() != other->nRows())
        throw new Exception("BlockVector: inconsistent dimensions in add!");

    for (unsigned int i = 0; i < nRows(); i++)
        block(i)->add(otherVector->block(i));

    updateNormInf();
}

void
BlockVector::
multiplyByScalar(const double& coeff)
{
    checkClosed();

    for (unsigned int i = 0; i < nRows(); i++)
        block(i)->multiplyByScalar(coeff);

    updateNormInf();
}

void
BlockVector::
dump(std::string namefile) const
{
    throw new Exception("Dump is not implemented for Block");
}

void
BlockVector::
softCopy(std::shared_ptr<aVector> other)
{
    if (other)
    {
        checkType(other, BLOCK);
        auto otherVector = dynamic_cast<BlockVector*>(other.get());

        resize(other->nRows());

        for (unsigned int i = 0; i < nRows(); i++)
            M_vectorGrid(i,0)->softCopy(otherVector->block(i));

        M_isOpen = otherVector->isOpen();
        M_normInf = otherVector->normInf();
    }
}

void
BlockVector::
hardCopy(std::shared_ptr<aVector> other)
{
    if (other)
    {
        checkType(other, BLOCK);
        auto otherVector = dynamic_cast<BlockVector*>(other.get());

        resize(other->nRows());

        for (unsigned int i = 0; i < nRows(); i++)
            M_vectorGrid(i,0)->hardCopy(otherVector->block(i));

        M_isOpen = otherVector->isOpen();
        M_normInf = otherVector->normInf();
    }
}

aVector*
BlockVector::
clone() const
{
    return new BlockVector(*this);
}

bool
BlockVector::
isZero() const
{
    for (unsigned int i = 0; i < nRows(); i++)
    {
        if (block(i) && !block(i)->isZero())
            return false;
    }
    return M_normInf < ZEROTHRESHOLD;
}

std::string
BlockVector::
getString(const char& delimiter) const
{
    throw new Exception("GetString not implement for BlockVector!");
}

std::shared_ptr<aVector>
BlockVector::
block(const unsigned int& iblock) const
{
    return M_vectorGrid(iblock,1);
}

std::shared_ptr<BlockVector>
BlockVector::
getSubvector(const unsigned int& ibegin, const unsigned int& iend) const
{
    std::shared_ptr<BlockVector> retVector(new BlockVector());

    unsigned int nrows = iend-ibegin+1;
    retVector->resize(iend-ibegin+1);

    for (unsigned int i = ibegin; i <= iend; i++)
        retVector->setBlock(i-ibegin,block(i));

    return retVector;
}

void
BlockVector::
setBlock(const unsigned int& iblock, std::shared_ptr<aVector> vector)
{
    checkOpen();
    M_vectorGrid(iblock,0) = vector;
    updateNormInf();
}

// BlockVector
// BlockVector::
// operator*(const double& coeff) const
// {
//     BlockVector ret;
//     ret.hardCopy(*this);
//
//     ret *= coeff;
//
//     return ret;
// }

double
BlockVector::
norm2() const
{
    if (nRows() == 0)
        return 0;

    double ret = 0;

    for (unsigned int i = 0; i < nRows(); i++)
    {
        double curNorm = block(i)->norm2();
        ret += curNorm * curNorm;
    }

    return sqrt(ret);
}

void
BlockVector::
updateNormInf()
{
    for (unsigned int i = 0; i < nRows(); i++)
    {
        if (block(i))
            M_normInf = M_normInf < block(i)->normInf() ? block(i)->normInf() : M_normInf;
    }
}

// BlockVector&
// BlockVector::
// operator*=(const double& coeff)
// {
//     for (unsigned int i = 0; i < M_nRows; i++)
//         block(i) *= coeff;
//
//     return *this;
// }
//
// BlockVector
// BlockVector::
// operator+(const BlockVector<InVectorType>& other) const
// {
//     BlockVector<InVectorType> ret;
//
//     if (M_nRows == 0)
//     {
//         ret.hardCopy(other);
//         return ret;
//     }
//
//     if (other.nRows() == 0)
//     {
//         ret.hardCopy(*this);
//         return ret;
//     }
//
//     if (M_nRows != other.M_nRows)
//     {
//         throw new Exception("Dimension of vectors being added is not consistent!");
//     }
//
//     ret.hardCopy(*this);
//
//     ret += other;
//     return ret;
// }
//
// BlockVector
// BlockVector::
// operator-(const BlockVector<InVectorType>& other) const
// {
//     BlockVector<InVectorType> ret;
//
//     if (M_nRows == 0)
//     {
//         ret.hardCopy(other);
//         return ret;
//     }
//
//     if (other.nRows() == 0)
//     {
//         ret.hardCopy(*this);
//         return ret;
//     }
//
//     if (M_nRows != other.M_nRows)
//     {
//         throw new Exception("Dimension of vectors being subtracted is not consistent!");
//     }
//
//     ret.hardCopy(*this);
//
//     ret -= other;
//     return ret;
// }
//
// BlockVector&
// BlockVector::
// operator+=(const BlockVector<InVectorType>& other)
// {
//     if (M_nRows == 0)
//     {
//         hardCopy(other);
//         return *this;
//     }
//
//     if (other.nRows() == 0)
//         return *this;
//
//     if (M_nRows != other.M_nRows)
//     {
//         throw new Exception("Dimension of vectors being added is not consistent!");
//     }
//
//     for (unsigned int i = 0; i < M_nRows; i++)
//         block(i) += other.block(i);
//
//     return *this;
// }
//
// BlockVector&
// BlockVector::
// operator-=(const BlockVector<InVectorType>& other)
// {
//     if (M_nRows == 0)
//     {
//         hardCopy(other);
//         (*this) *= (-1);
//         return *this;
//     }
//
//     if (other.nRows() == 0)
//         return *this;
//
//     if (M_nRows != other.M_nRows)
//     {
//         throw new Exception("Dimension of vectors being subtracted is not consistent!");
//     }
//
//     for (unsigned int i = 0; i < M_nRows; i++)
//         block(i) -= other.block(i);
//
//     return *this;
// }

// void
// BlockVector::
// hardCopy(const BlockVector& other)
// {
//     // resize(other.M_nRows);
//     //
//     // for (unsigned int i = 0; i < M_nRows; i++)
//     //     block(i).hardCopy(other.block(i));
// }
//
// void
// BlockVector::
// softCopy(const BlockVector& other)
// {
//     // resize(other.M_nRows);
//     // for (unsigned int i = 0; i < M_nRows; i++)
//     // {
//     //     block(i).softCopy(other.block(i));
//     // }
// }

}
