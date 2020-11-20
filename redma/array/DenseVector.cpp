#include "DenseVector.hpp"

namespace RedMA
{

DenseVector::
DenseVector() :
  M_vector(nullptr)
{

}

DenseVector::
DenseVector(const DenseVector& vector)
{
    M_vector.reset(new DENSEVECTOR(*vector.M_vector));
}

void
DenseVector::
add(std::shared_ptr<aVector> other)
{
    if (isZero())
    {
        deepCopy(other);
        return;
    }

    if (!other->isZero())
    {
        shp<DenseVector> otherVector = convert<DenseVector>(other);
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
shallowCopy(std::shared_ptr<aDataWrapper> other)
{
    if (other)
    {
        auto otherVector = convert<DenseVector>(other);
        setVector(otherVector->M_vector);
    }
}

void
DenseVector::
deepCopy(std::shared_ptr<aDataWrapper> other)
{
    if (other)
    {
        auto otherVector = convert<DenseVector>(other);
        std::shared_ptr<DENSEVECTOR> newVector
            (new DENSEVECTOR(*otherVector->M_vector));
        setVector(newVector);
    }
}

DenseVector*
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
isZero() const
{
    return M_vector == nullptr;
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
toVectorEpetraPtr(std::shared_ptr<Epetra_Comm> comm) const
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
    }
}


};
