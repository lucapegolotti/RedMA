#include "DenseVector.hpp"

namespace RedMA
{

DenseVector::
DenseVector() :
  M_vector(nullptr)
{

}

void
DenseVector::
add(shp<aVector> other)
{
    if (isZero())
    {
        deepCopy(other);
        return;
    }

    if (!other->isZero())
    {
        if (other->nRows() != nRows())
            throw new Exception("[DenseVector::add] inconsistent dimensions of vectors");

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
shallowCopy(shp<aDataWrapper> other)
{
    if (other)
    {
        auto otherVector = convert<DenseVector>(other);
        setVector(otherVector->M_vector);
    }
}

void
DenseVector::
deepCopy(shp<aDataWrapper> other)
{
    if (other)
    {
        auto otherVector = convert<DenseVector>(other);
        shp<DENSEVECTOR> newVector
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
        shp<DENSEVECTOR> newVector
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
setData(shp<void> data)
{
    setVector(spcast<DENSEVECTOR>(data));
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

shp<LifeV::VectorEpetra>
DenseVector::
toVectorEpetraPtr(shp<Epetra_Comm> comm) const
{
    // note: we dont care about parallelism because we are assume that we are serial
    using namespace LifeV;
    unsigned int N = M_vector->Length();

    shp<MapEpetra> rangeMap;
    rangeMap.reset(new MapEpetra(N, N, 0, comm));

    shp<VectorEpetra> retVec(new VectorEpetra(*rangeMap));

    for (unsigned int i = 0; i < N; i++)
        (*retVec)[i] = (*M_vector)(i);

    return retVec;
}

void
DenseVector::
setVector(shp<DENSEVECTOR> vector)
{
    if (vector)
    {
        M_vector = vector;
        this->M_nRows = M_vector->Length();
    }
}

shp<DENSEVECTOR>
DenseVector::
getVector()
{
    return M_vector;
}

shp<void>
DenseVector::
data() const
{
    return M_vector;
}

};
