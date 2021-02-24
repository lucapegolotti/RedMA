#include "DistributedVector.hpp"

namespace RedMA
{

DistributedVector::
DistributedVector() :
  M_vector(nullptr)
{

}

DistributedVector::
DistributedVector(const DistributedVector& vector)
{
    M_vector.reset(new VECTOREPETRA(*vector.M_vector));
}

void
DistributedVector::
add(shp<aVector> other)
{
    if (other->isZero())
        return;

    if (other->type() == DENSE)
        other = convertDenseVector(spcast<DenseVector>(other),commPtr());

    auto otherVector = convert<DistributedVector>(other);

    if (isZero())
    {
        deepCopy(other);
        return;
    }

    if (!other->isZero())
    {
        if (other->nRows() != nRows())
            throw new Exception("[DistributedVector::add] inconsistent dimensions of vectors");
        (*M_vector) += *otherVector->getVector();
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
    M_vector->spy(namefile);
}

void
DistributedVector::
shallowCopy(shp<aDataWrapper> other)
{
    if (other)
    {
        auto otherVector = convert<DistributedVector>(other);
        setVector(otherVector->M_vector);
    }
}

void
DistributedVector::
deepCopy(shp<aDataWrapper> other)
{
    if (other)
    {
        auto otherVector = convert<DistributedVector>(other);
        shp<VECTOREPETRA> newVector;
        if (otherVector->M_vector)
            newVector.reset(new VECTOREPETRA(*otherVector->M_vector));
        setVector(newVector);
    }
}

DistributedVector*
DistributedVector::
clone() const
{
    DistributedVector* retVector = new DistributedVector();
    if (M_vector)
    {
        shp<VECTOREPETRA> newVector
            (new VECTOREPETRA(*M_vector));
        retVector->setVector(newVector);
    }
    return retVector;
}

bool
DistributedVector::
isZero() const
{
    return M_vector == nullptr;
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

void
DistributedVector::
setData(shp<void> data)
{
    setVector(spcast<VECTOREPETRA>(data));
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
    return *toDenseVectorPtr();
}

shp<DenseVector>
DistributedVector::
toDenseVectorPtr() const
{
    shp<DenseVector> retVec(new DenseVector());

    if (M_vector)
    {
        auto mapPtr = M_vector->mapPtr();
        shp<DENSEVECTOR> innerVector(new DENSEVECTOR(mapPtr->mapSize()));

        for (unsigned int i = 0; i < innerVector->Length(); i++)
            (*innerVector)(i) = M_vector->operator[](i);

        retVec->setVector(innerVector);
    }

    return retVec;
}

shp<DistributedVector>
DistributedVector::
convertDenseVector(shp<DenseVector> denseVector, shp<Epetra_Comm> comm)
{
    using namespace LifeV;
    shp<DistributedVector> retVec(new DistributedVector());

    if (comm->MyPID() != 0)
        throw new Exception("convertDenseVector does not support more than one proc");

    unsigned int length = denseVector->nRows();
    shp<MapEpetra> map(new MapEpetra(length, length, 0, comm));

    shp<VECTOREPETRA> innerVec(new VECTOREPETRA(*map));
    innerVec->zero();

    for (unsigned int i = 0; i < length; i++)
    {
        double value = spcast<DENSEVECTOR>(denseVector->data())->operator()(i);
        innerVec->operator[](i) = value;
    }
    shp<DistributedVector> retVector(new DistributedVector());
    retVector->setVector(innerVec);

    return retVector;
}


void
DistributedVector::
setVector(shp<VECTOREPETRA> vector)
{
    if (vector)
    {
        auto mapPtr = vector->mapPtr();

        M_vector = vector;
        this->M_nRows = mapPtr->mapSize();
    }
}

shp<void>
DistributedVector::
data() const
{
    return M_vector;
}

};
