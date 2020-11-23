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
add(std::shared_ptr<aVector> other)
{
    if (other->isZero())
        return;

    if (other->type() == DENSE)
        other = convertDenseVector(std::static_pointer_cast<DenseVector>(other),commPtr());

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
    throw new Exception("Dump not implemented for DistributedVector");
}

void
DistributedVector::
shallowCopy(std::shared_ptr<aDataWrapper> other)
{
    if (other)
    {
        auto otherVector = convert<DistributedVector>(other);
        setVector(otherVector->M_vector);
    }
}

void
DistributedVector::
deepCopy(std::shared_ptr<aDataWrapper> other)
{
    if (other)
    {
        auto otherVector = convert<DistributedVector>(other);
        std::shared_ptr<VECTOREPETRA> newVector;
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
        std::shared_ptr<VECTOREPETRA> newVector
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

std::shared_ptr<void>
DistributedVector::
data() const
{
    return M_vector;
}

void
DistributedVector::
setData(std::shared_ptr<void> data)
{
    M_vector = std::static_pointer_cast<VECTOREPETRA>(data);
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

std::shared_ptr<DenseVector>
DistributedVector::
toDenseVectorPtr() const
{
    std::shared_ptr<DenseVector> retVec(new DenseVector());

    if (M_vector)
    {
        auto mapPtr = M_vector->mapPtr();
        std::shared_ptr<DENSEVECTOR> innerVector(new DENSEVECTOR(mapPtr->mapSize()));

        for (unsigned int i = 0; i < innerVector->Length(); i++)
            (*innerVector)(i) = M_vector->operator[](i);

        retVec->setVector(innerVector);
    }

    return retVec;
}

// std::shared_ptr<DistributedVector>
// DistributedVector::
// convertDenseVector(std::shared_ptr<DenseVector> denseVector, std::shared_ptr<Epetra_Comm> comm)
// {
//     using namespace LifeV;
//     std::shared_ptr<DistributedVector> retVec(new DistributedVector());
//
//     if (comm->MyPID() != 0)
//         throw new Exception("convertDenseVector does not support more than one proc");
//
//     unsigned int length = denseVector->nRows();
//     std::cout << "length = " << length << std::endl << std::flush;
//     std::shared_ptr<MapEpetra> map(new MapEpetra(length, length, 0, comm));
//
//     retVec->data().reset(new VECTOREPETRA(*map));
//
//     for (unsigned int i = 0; i < length; i++)
//     {
//         static_cast<VECTOREPETRA*>(retVec->data().get())->operator[](i) =
//         static_cast<VECTOREPETRA*>(denseVector->data().get())->operator()(i);
//     }
//
//     return retVec;
// }

std::shared_ptr<DistributedVector>
DistributedVector::
convertDenseVector(std::shared_ptr<DenseVector> denseVector, std::shared_ptr<Epetra_Comm> comm)
{
    using namespace LifeV;
    std::shared_ptr<DistributedVector> retVec(new DistributedVector());

    if (comm->MyPID() != 0)
        throw new Exception("convertDenseVector does not support more than one proc");

    unsigned int length = denseVector->nRows();
    std::shared_ptr<MapEpetra> map(new MapEpetra(length, length, 0, comm));

    std::shared_ptr<VECTOREPETRA> innerVec(new VECTOREPETRA(*map));
    innerVec->zero();

    for (unsigned int i = 0; i < length; i++)
    {
        double value = std::static_pointer_cast<DENSEVECTOR>(denseVector->data())->operator()(i);
        innerVec->operator[](i) = value;
    }
    std::shared_ptr<DistributedVector> retVector(new DistributedVector());
    retVector->setVector(innerVec);

    return retVector;
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
    }
}

std::shared_ptr<DistributedVector>
epetraToDistributed(std::shared_ptr<VECTOREPETRA> vector)
{
    std::shared_ptr<DistributedVector> retVec(new DistributedVector());
    retVec->setData(vector);
    return retVec;
}

};
