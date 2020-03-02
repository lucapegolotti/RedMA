#include "VectorEp.hpp"

namespace RedMA
{

VectorEp::
VectorEp() :
  M_vector(nullptr)
{

}

VectorEp
VectorEp::
operator+(const VectorEp& other)
{
    VectorEp vec;

    if (!M_vector)
    {
        hardCopy(other);
        return other;
    }

    vec.hardCopy(*this);
    vec += other;
    return vec;
}

VectorEp
VectorEp::
operator-(const VectorEp& other)
{
    VectorEp vec;

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

VectorEp&
VectorEp::
operator+=(const VectorEp& other)
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

VectorEp&
VectorEp::
operator-=(const VectorEp& other)
{
    if (!M_vector)
    {
        hardCopy(other);
        *M_vector *= (-1.0);
        return *this;
    }

    if (other.data())
        *M_vector -= *other.data();
    return *this;
}

VectorEp&
VectorEp::
operator=(const std::shared_ptr<VECTOREPETRA>& other)
{
    M_vector = other;
    return *this;
}

VectorEp&
VectorEp::
operator*=(const double& coeff)
{
    if (M_vector)
        *M_vector *= coeff;

    return *this;
}

void
VectorEp::
hardCopy(const VectorEp& other)
{
    if (other.data())
        M_vector.reset(new VECTOREPETRA(*other.data()));
}

void
VectorEp::
softCopy(const VectorEp& other)
{
    M_vector = other.data();
}

double
VectorEp::
norm2() const
{
    double mynorm = 0;
    if (M_vector)
        mynorm += M_vector->norm2();

    return mynorm;
}

std::shared_ptr<VECTOREPETRA>&
VectorEp::
data()
{
    return M_vector;
}

std::shared_ptr<VECTOREPETRA>
VectorEp::
data() const
{
    return M_vector;
}

std::string
VectorEp::
getString(const char& delimiter) const
{
    std::string ret = "";
    // we reduce the vector to processor 0
    VECTOREPETRA redVec(*M_vector, 0);

    if (redVec.epetraVector().Comm().MyPID())
        return ret;

    const double* values = redVec.epetraVector()[0];
    for (unsigned int i = 0; i < redVec.epetraVector().GlobalLength(); ++i)
    {
        ret += std::to_string(values[i]);
        if (i != redVec.epetraVector().GlobalLength()-1)
            ret += delimiter;
    }
    return ret;
}

// here we assume that the vector is stacked 1st component, 2nd component ecc
double
VectorEp::
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


};
