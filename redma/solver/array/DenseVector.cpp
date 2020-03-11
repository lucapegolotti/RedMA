#include "DenseVector.hpp"

namespace RedMA
{

DenseVector::
DenseVector() :
  M_vector(nullptr)
{

}

DenseVector
DenseVector::
operator+(const DenseVector& other)
{
    DenseVector vec;

    if (!M_vector)
    {
        hardCopy(other);
        return other;
    }

    vec.hardCopy(*this);
    vec += other;
    return vec;
}

DenseVector
DenseVector::
operator-(const DenseVector& other)
{
    DenseVector vec;

    if (!M_vector)
    {
        hardCopy(other);
        *this *= -1.0;
        return other;
    }

    vec.hardCopy(*this);
    vec -= other;
    return vec;
}

DenseVector&
DenseVector::
operator+=(const DenseVector& other)
{
    if (!M_vector)
    {
        hardCopy(other);
        return *this;
    }

    if (other.data())
        *M_vector += *other.data();
    return *this;
}

DenseVector&
DenseVector::
operator-=(const DenseVector& other)
{
    if (!M_vector)
    {
        hardCopy(other);
        M_vector->Scale(-1.0);
        return *this;
    }

    if (other.data())
    {
        for (unsigned int i = 0; i < M_vector->Length(); i++)
            (*M_vector)[i] -= (*other.data())[i];
    }
    return *this;
}

DenseVector&
DenseVector::
operator=(const std::shared_ptr<DENSEVECTOR>& other)
{
    M_vector = other;
}

DenseVector&
DenseVector::
operator*=(const double& coeff)
{
    if (M_vector)
        M_vector->Scale(coeff);

    return *this;
}

void
DenseVector::
hardCopy(const DenseVector& other)
{
    if (other.data())
        M_vector.reset(new DENSEVECTOR(*other.data()));
}

void
DenseVector::
softCopy(const DenseVector& other)
{
    M_vector = other.data();
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

std::shared_ptr<DENSEVECTOR>&
DenseVector::
data()
{
    return M_vector;
}

std::shared_ptr<DENSEVECTOR>
DenseVector::
data() const
{
    return M_vector;
}

std::string
DenseVector::
getString(const char& delimiter) const
{
    std::string ret = "";

    for (unsigned int i = 0; i < M_vector->Length(); ++i)
    {
        ret += std::to_string((*M_vector)[i]);
        if (i != M_vector->Length()-1)
            ret += delimiter;
    }
    return ret;
}

unsigned int
DenseVector::
getNumRows() const
{
    if (!M_vector)
        throw new Exception("DenseVector::getNumRows() inner vector not set");

    return M_vector->Length();
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

};
