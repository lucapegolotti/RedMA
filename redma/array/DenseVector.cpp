#include "DenseVector.hpp"

namespace RedMA
{

DenseVector::
DenseVector() :
  aVector(DENSE),
  M_vector(nullptr)
{

}

DenseVector::
DenseVector(const DenseVector& vector) :
  aVector(DENSE)
{
    M_vector.reset(new DENSEVECTOR(*vector.M_vector));
}

void
DenseVector::
add(std::shared_ptr<aVector> other)
{
    checkType(other, DENSE);

    if (isZero())
    {
        hardCopy(other);
        return;
    }

    if (!other->isZero())
    {
        DenseVector* otherVector = dynamic_cast<DenseVector*>(other.get());
        (*M_vector) += (*static_cast<DENSEVECTOR*>(otherVector->data().get()));
    }
}

void
DenseVector::
multiplyByScalar(const double& coeff)
{
    if (!isZero())
        M_vector->Scale(coeff);
}

void
DenseVector::
softCopy(std::shared_ptr<aVector> other)
{
    if (other)
    {
        checkType(other, DENSE);
        auto otherMatrix = dynamic_cast<DenseVector*>(other.get());
        setVector(otherMatrix->M_vector);
    }
}

void
DenseVector::
hardCopy(std::shared_ptr<aVector> other)
{
    if (other)
    {
        checkType(other, DENSE);
        auto otherVector = dynamic_cast<DenseVector*>(other.get());
        std::shared_ptr<DENSEVECTOR> newVector
            (new DENSEVECTOR(*otherVector->M_vector));
        setVector(newVector);
    }
}

aVector*
DenseVector::
clone() const
{
    DenseVector* retVector = new DenseVector();
    if (M_vector)
    {
        std::shared_ptr<DENSEVECTOR> newVector
            (new DENSEVECTOR(*M_vector));
        retVector->setVector(newVector);
    }
    return retVector;
}

bool
DenseVector::
isZero()
{
    if (!M_vector)
        return true;
    // we recompute just to make sure
    if (M_normInf < ZEROTHRESHOLD)
        M_normInf = M_vector->NormInf();
    return M_normInf < ZEROTHRESHOLD;
}

double
DenseVector::
norm2() const
{
    double mynorm = 0;
    if (M_vector)
        mynorm += M_vector->Norm2();

    return mynorm;
}

void
DenseVector::
setData(std::shared_ptr<void> data)
{
    M_vector = std::static_pointer_cast<DENSEVECTOR>(data);
}

std::shared_ptr<void>
DenseVector::
data() const
{
    return M_vector;
}

std::string
DenseVector::
getString(const char& delimiter) const
{
    std::ostringstream streamObj;
    // streamObj << std::scientific;
    streamObj << std::setprecision(16);
    streamObj << "";

    if (M_vector)
    {
        for (unsigned int i = 0; i < M_vector->Length(); ++i)
        {
            if (std::abs((*M_vector)[i]) > ZEROTHRESHOLD)
                streamObj << (*M_vector)[i];
            else
                streamObj << 0.0;
            if (i != M_vector->Length()-1)
                streamObj << delimiter;
        }
    }
    return streamObj.str();
}

void
DenseVector::
dump(std::string filename) const
{
    std::ofstream outfile(filename);

    if (M_vector)
    {
        unsigned int N = M_vector->Length();

        for (unsigned int i = 0; i < N; i++)
        {
            outfile << (*M_vector)(i);
            outfile << "\n";
        }
    }

    outfile.close();
}

std::shared_ptr<LifeV::VectorEpetra>
DenseVector::
toVectorEpetra(std::shared_ptr<Epetra_Comm> comm) const
{
    // note: we dont care about parallelism because we are assume that we are serial
    using namespace LifeV;
    unsigned int N = M_vector->Length();

    std::shared_ptr<MapEpetra> rangeMap;
    rangeMap.reset(new MapEpetra(N, N, 0, comm));

    std::shared_ptr<VectorEpetra> retVec(new VectorEpetra(*rangeMap));

    for (unsigned int i = 0; i < N; i++)
        (*retVec)[i] = (*M_vector)(i);

    return retVec;
}

void
DenseVector::
setVector(std::shared_ptr<DENSEVECTOR> vector)
{
    if (vector)
    {
        M_vector = vector;
        this->M_nRows = M_vector->Length();
        this->M_normInf = M_vector->NormInf();
    }
}


};
