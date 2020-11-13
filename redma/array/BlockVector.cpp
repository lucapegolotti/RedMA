#include "BlockVector.hpp"

namespace RedMA
{

BlockVector::
BlockVector() :
  aVector(BLOCK)
  // M_isOpen(true)
{
}

BlockVector::
BlockVector(const BlockVector& other) :
  aVector(BLOCK)
  // M_isOpen(other.M_isOpen)
{
    // other.checkClosed();
    resize(other.nRows());

    for (unsigned int i = 0; i < other.nRows(); i++)
    {
        setBlock(i, other.block(i));
    }
}


BlockVector::
BlockVector(const unsigned int& nRows) :
  aVector(BLOCK)
  // M_isOpen(true)
{
    resize(nRows);
}

void
BlockVector::
findComm()
{
    M_comm = nullptr;

    for (unsigned int i = 0; i < nRows(); i++)
    {
        if (block(i) && M_comm == nullptr)
        {
            if (block(i)->type() == BLOCK)
                M_comm = std::static_pointer_cast<BlockVector>(block(i))->commPtr();
            else if (block(i)->type() == DISTRIBUTED)
            {
                M_comm = std::static_pointer_cast<DistributedVector>(block(i))->commPtr();
            }
        }
    }
}

bool
BlockVector::
globalTypeIs(Datatype type)
{
    for (unsigned int i = 0; i < nRows(); i++)
    {
        if (block(i))
        {
            if (!block(i)->isZero() && block(i)->type() == BLOCK)
            {
                if (!std::static_pointer_cast<BlockVector>(block(i))->globalTypeIs(type))
                    return false;
            }
            else if (!block(i)->isZero() && block(i)->type() != type)
            {
                return false;
            }
        }
    }
    return true;
}

void
BlockVector::
copyPattern(std::shared_ptr<BlockVector> other, bool verbose)
{
    throw new Exception("This function does not work");
    for (unsigned int i = 0; i < nRows(); i++)
    {
        if (block(i)->type() == BLOCK)
        {
            std::static_pointer_cast<BlockVector>(block(i))->copyPattern(std::static_pointer_cast<BlockVector>(other->block(i)),false);
        }
        else if (block(i)->type() != other->block(i)->type())
        {
            if (other->block(i)->type() == DENSE && block(i)->type() == DISTRIBUTED)
            {
                setBlock(i,std::static_pointer_cast<DistributedVector>(block(i))->toDenseVectorPtr());
            }
            else
                throw new Exception("copyPattern: case not implemented");
        }
    }
}

std::shared_ptr<BlockVector>
BlockVector::
convertInnerTo(Datatype type, std::shared_ptr<Epetra_Comm> comm)
{
    if (type != DISTRIBUTED)
        throw new Exception("convertInnerTo implemented only for distributed output");

    if (comm == nullptr)
    {
        findComm();
        comm = M_comm;
    }

    if (comm == nullptr)
        throw new Exception("No communicator in block vector!");

    std::shared_ptr<BlockVector> retVec(new BlockVector(nRows()));

    for (unsigned int i = 0; i < nRows(); i++)
    {
        retVec->setBlock(i,block(i));
        if (block(i))
        {
            if (!block(i)->isZero() && block(i)->type() == BLOCK)
            {
                retVec->setBlock(i,std::static_pointer_cast<BlockVector>(block(i))->convertInnerTo(type,comm));
            }
            else if (!block(i)->isZero() && block(i)->type() != type)
            {
                if (block(i)->type() != DENSE)
                    throw new Exception("All inner vectors which are not sparse must be dense!");
                retVec->setBlock(i,DistributedVector::convertDenseVector(std::static_pointer_cast<DenseVector>(block(i)),comm));
            }
        }
    }

    return retVec;
}

void
BlockVector::
resize(const unsigned int& nRows)
{
    M_nRows = nRows;
    M_vectorGrid.resize(nRows,1);

    for (unsigned int i = 0; i < nRows; i++)
        M_vectorGrid(i,0).reset(new DenseVector());
}

void
BlockVector::
add(std::shared_ptr<aVector> other)
{
    if (other->nRows() == 0)
        return;

    if (nRows() == 0)
    {
        hardCopy(other);
        return;
    }
    std::shared_ptr<BlockVector> otherVector = std::static_pointer_cast<BlockVector>(other);

    if (nRows() != other->nRows())
        throw new Exception("BlockVector: inconsistent dimensions in add!");

    for (unsigned int i = 0; i < nRows(); i++)
    {
        if (!block(i)->isZero())
        {
            if (otherVector && !otherVector->block(i)->isZero())
            {
                block(i)->add(otherVector->block(i));
            }
        }
        else
            if (otherVector->block(i))
                setBlock(i,std::shared_ptr<aVector>(otherVector->block(i)->cloneVector()));
    }

    updateNormInf();
}

void
BlockVector::
multiplyByScalar(const double& coeff)
{
    for (unsigned int i = 0; i < nRows(); i++)
    if (block(i))
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

        // M_isOpen = otherVector->isOpen();
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
        std::shared_ptr<BlockVector> otherVector = std::static_pointer_cast<BlockVector>(other);

        resize(other->nRows());

        for (unsigned int i = 0; i < nRows(); i++)
        {
            if (otherVector->block(i))
                M_vectorGrid(i,0).reset(otherVector->block(i)->cloneVector());
        }

        // M_isOpen = otherVector->isOpen();
        M_normInf = otherVector->normInf();
    }
}

aVector*
BlockVector::
cloneVector() const
{
    BlockVector* retVector = new BlockVector(nRows());
    for (unsigned int i = 0; i < nRows(); i++)
        retVector->setBlock(i,std::shared_ptr<aVector>(block(i)->cloneVector()));
    return retVector;
}

bool
BlockVector::
isZero()
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
    return M_vectorGrid(iblock,0);
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
