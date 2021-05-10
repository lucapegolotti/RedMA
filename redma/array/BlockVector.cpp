#include "BlockVector.hpp"

namespace RedMA
{

BlockVector::
BlockVector()
{
}

BlockVector::
BlockVector(const unsigned int& nRows)
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
                M_comm = convert<BlockVector>(block(i))->commPtr();
            else if (block(i)->type() == DISTRIBUTED)
            {
                M_comm = convert<DistributedVector>(block(i))->commPtr();
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
                if (!convert<BlockVector>(block(i))->globalTypeIs(type))
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

shp<BlockVector>
BlockVector::
convertInnerTo(Datatype type, shp<Epetra_Comm> comm)
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

    shp<BlockVector> retVec(new BlockVector(nRows()));

    for (unsigned int i = 0; i < nRows(); i++)
    {
        retVec->setBlock(i,block(i));
        if (block(i))
        {
            if (!block(i)->isZero() && block(i)->type() == BLOCK)
            {
                retVec->setBlock(i,spcast<BlockVector>(block(i))->convertInnerTo(type,comm));
            }
            else if (!block(i)->isZero() && block(i)->type() != type)
            {
                if (block(i)->type() != DENSE)
                    throw new Exception("All inner vectors which are not sparse must be dense!");
                retVec->setBlock(i,DistributedVector::convertDenseVector(spcast<DenseVector>(block(i)),comm));
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
add(shp<aVector> other)
{
    if (other->nRows() == 0)
        return;

    if (nRows() == 0)
    {
        deepCopy(other);
        return;
    }
    shp<BlockVector> otherVector = convert<BlockVector>(other);

    if (nRows() != other->nRows())
        throw new Exception("BlockVector: inconsistent dimensions in add!");

    for (unsigned int i = 0; i < nRows(); i++)
    {
        if (!block(i)->isZero())
        {
            if (otherVector && !otherVector->block(i)->isZero())
            {
                if (otherVector->block(i)->type() == DISTRIBUTED &&
                    block(i)->type() == DENSE)
                {
                    shp<DistributedVector> asDistributed = convert<DistributedVector>(otherVector->block(i));
                    shp<DenseVector> asDense = asDistributed->toDenseVectorPtr();
                    block(i)->add(asDense);
                }
                else
                    block(i)->add(otherVector->block(i));
            }
        }
        else
            if (otherVector->block(i))
            {
                shp<aDataWrapper> blck(otherVector->block(i)->clone());
                setBlock(i,spcast<aVector>(blck));
            }
    }
}

void
BlockVector::
multiplyByScalar(const double& coeff)
{
    for (unsigned int i = 0; i < nRows(); i++)
    if (block(i))
        block(i)->multiplyByScalar(coeff);
}

void
BlockVector::
dump(std::string namefile) const
{
    throw new Exception("Dump is not implemented for Block");
}

void
BlockVector::
shallowCopy(shp<aDataWrapper> other)
{
    if (other && !other->isZero())
    {
        auto otherVector = convert<BlockVector>(other);

        resize(otherVector->nRows());

        for (unsigned int i = 0; i < nRows(); i++)
            M_vectorGrid(i,0) = otherVector->block(i);

    }
}

void
BlockVector::
deepCopy(shp<aDataWrapper> other)
{
    if (other && !other->isZero())
    {
        shp<BlockVector> otherVector = spcast<BlockVector>(other);

        resize(otherVector->nRows());

        for (unsigned int i = 0; i < nRows(); i++)
        {
            if (otherVector->block(i))
                M_vectorGrid(i,0).reset(static_cast<aVector*>(otherVector->block(i)->clone()));
        }
    }
}

BlockVector*
BlockVector::
clone() const
{
    BlockVector* retVector = new BlockVector(nRows());
    for (unsigned int i = 0; i < nRows(); i++)
        retVector->setBlock(i,shp<aVector>(static_cast<aVector*>(block(i)->clone())));
    return retVector;
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
    return true;
}

std::string
BlockVector::
getString(const char& delimiter) const
{
    throw new Exception("GetString not implement for BlockVector!");
}

shp<aVector>
BlockVector::
block(const unsigned int& iblock) const
{
    if (iblock >= nRows())
        throw new Exception("BlockVector: iblock > nRows()!");
    return M_vectorGrid(iblock,0);
}

shp<BlockVector>
BlockVector::
getSubvector(const unsigned int& ibegin, const unsigned int& iend) const
{
    shp<BlockVector> retVector(new BlockVector());

    unsigned int nrows = iend-ibegin+1;
    retVector->resize(iend-ibegin+1);

    for (unsigned int i = ibegin; i <= iend; i++)
        retVector->setBlock(i-ibegin,block(i));

    return retVector;
}

void
BlockVector::
setBlock(const unsigned int& iblock, shp<aVector> vector)
{
    if (iblock >= nRows())
        throw new Exception("[BlockVector::setBlock] iblock >= nRows()!!");
    M_vectorGrid(iblock,0) = vector;
}

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

unsigned int
BlockVector::
level()
{
    for (unsigned int i = 0; i < nRows(); i++)
    {
        if (block(i)->type() == BLOCK)
            return convert<BlockVector>(block(i))->level() + 1;
    }
    return 1;
}

}
