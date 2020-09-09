#include "DistributedVector.hpp"

namespace RedMA
{

DistributedVector::
DistributedVector() :
  aVector(DISTRIBUTED),
  M_vector(nullptr)
{

}

void
DistributedVector::
add(std::shared_ptr<aVector> other)
{
    checkType(other, DISTRIBUTED);

    if (isZero())
    {
        hardCopy(other);
        return;
    }

    if (!other->isZero())
    {
        DistributedVector* otherVector = dynamic_cast<DistributedVector*>(other.get());
        (*M_vector) += (*static_cast<VECTOREPETRA*>(otherVector->data().get()));
    }
}

void
DistributedVector::
multiplyByScalar(const double& coeff)
{
    if (!isZero())
        (*M_vector) *= coeff;
}

void
DistributedVector::
dump(std::string namefile) const
{
    throw new Exception("Dump not implemented for DistributedVector");
}

void
DistributedVector::
softCopy(std::shared_ptr<aVector> other)
{
    if (other)
    {
        checkType(other, DISTRIBUTED);
        auto otherVector = dynamic_cast<DistributedVector*>(other.get());
        setVector(otherVector->M_vector);
    }
}

void
DistributedVector::
hardCopy(std::shared_ptr<aVector> other)
{
    if (other)
    {
        checkType(other, DISTRIBUTED);
        auto otherVector = dynamic_cast<DistributedVector*>(other.get());
        std::shared_ptr<VECTOREPETRA> newVector
            (new VECTOREPETRA(*otherVector->M_vector));
        setVector(newVector);
    }
}

aVector*
DistributedVector::
clone() const
{
    return new DistributedVector(*this);
}

bool
DistributedVector::
isZero() const
{
    if (!M_vector)
        return true;
    return M_normInf < ZEROTHRESHOLD;
}

double
DistributedVector::
norm2() const
{
    double mynorm = 0;
    if (M_vector)
        mynorm += M_vector->norm2();

    return mynorm;
}

std::shared_ptr<void>
DistributedVector::
data() const
{
    return M_vector;
}

std::string
DistributedVector::
getString(const char& delimiter) const
{
    std::ostringstream streamObj;
    // streamObj << std::scientific;
    streamObj << std::setprecision(16);
    streamObj << "";
    if (M_vector)
    {
        // we reduce the vector to processor 0
        VECTOREPETRA redVec(*M_vector, 0);

        if (redVec.epetraVector().Comm().MyPID())
            return streamObj.str();

        const double* values = redVec.epetraVector()[0];
        for (unsigned int i = 0; i < redVec.epetraVector().GlobalLength(); ++i)
        {
            if (std::abs(values[i]) > ZEROTHRESHOLD)
                streamObj << values[i];
            else
                streamObj << 0.;

            if (i != redVec.epetraVector().GlobalLength()-1)
                streamObj << delimiter;
        }
    }
    return streamObj.str();
}

// here we assume that the vector is stacked 1st component, 2nd component ecc
double
DistributedVector::
maxMagnitude3D() const
{
    double retval = -1;

    // we reduce the vector to processor 0
    VECTOREPETRA redVec(*M_vector, 0);

    if (redVec.epetraVector().GlobalLength() % 3 != 0)
        throw new Exception("norm2_3D: this is not a 3D vector!");

    unsigned int N = redVec.epetraVector().GlobalLength() / 3;

    if (redVec.epetraVector().Comm().MyPID())
        return retval;

    const double* values = redVec.epetraVector()[0];
    for (unsigned int i = 0; i < N; i++)
    {
        double curval = 0;
        for (unsigned int j = 0; j < 3; j++)
            curval += values[i + N*j] * values[i + N*j];
        curval = std::sqrt(curval);

        retval = retval < curval ? curval : retval;
    }

    return retval;
}

DenseVector
DistributedVector::
toDenseVector() const
{
    DenseVector retVec;

    if (M_vector)
    {
        auto mapPtr = M_vector->mapPtr();
        std::shared_ptr<DENSEVECTOR> innerVector(new DENSEVECTOR(mapPtr->mapSize()));

        for (unsigned int i = 0; i < innerVector->Length(); i++)
            (*innerVector)(i) = M_vector->operator[](i);

        retVec.data() = innerVector;
    }

    return retVec;
}

DistributedVector
DistributedVector::
convertDenseVector(DenseVector denseVector, std::shared_ptr<Epetra_Comm> comm)
{
    using namespace LifeV;
    DistributedVector retVec;

    if (comm->MyPID() != 0)
        throw new Exception("convertDenseVector does not support more than one proc");

    unsigned int length = denseVector.nRows();

    std::shared_ptr<MapEpetra> map(new MapEpetra(length, length, 0, comm));

    retVec.data().reset(new VECTOREPETRA(*map));

    for (unsigned int i = 0; i < length; i++)
    {
        static_cast<VECTOREPETRA*>(retVec.data().get())->operator[](i) =
        static_cast<VECTOREPETRA*>(denseVector.data().get())->operator()(i);
    }

    return retVec;
}

void
DistributedVector::
setVector(std::shared_ptr<VECTOREPETRA> vector)
{
    if (vector)
    {
        auto mapPtr = vector->mapPtr();

        M_vector = vector;
        this->M_nRows = mapPtr->mapSize();
        this->M_normInf = M_vector->normInf();
    }
}

};
