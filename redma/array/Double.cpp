#include "Double.hpp"

namespace RedMA
{

Double::
Double()
{
    M_double = 0;
}

Double
Double::
operator+(const Double& other)
{
    Double retDouble;
    retDouble.data() = M_double + other.data();
    return retDouble;
}

Double&
Double::
operator+=(const Double& other)
{
    M_double += other.data();
    return *this;
}

Double&
Double::
operator-=(const Double& other)
{
    M_double -= other.data();
    return *this;
}

Double&
Double::
operator*=(const double& coeff)
{
    M_double *= coeff;
    return *this;
}

Double&
Double::
operator=(const double& other)
{
    M_double = other;
    return *this;
}

Double
Double::
operator*(const Double& other)
{
    Double retDouble;
    retDouble.data() = M_double * other.data();
    return retDouble;
}


void
Double::
hardCopy(const Double& other)
{
    M_double = other.data();
}

void
Double::
softCopy(const Double& other)
{
    M_double = other.data();
}

double
Double::
norm2() const
{
    return M_double;
}

double&
Double::
data()
{
    return M_double;
}

double
Double::
data() const
{
    return M_double;
}

}
