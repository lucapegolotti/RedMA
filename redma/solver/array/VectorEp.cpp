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

VectorEp&
VectorEp::
operator+=(const VectorEp& other)
{
    if (!M_vector)
    {
        hardCopy(other);
        return *this;
    }

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

    *M_vector -= *other.data();
    return *this;
}

VectorEp&
VectorEp::
operator=(const std::shared_ptr<VECTOREPETRA>& other)
{
    std::cout << "Operator = could be buggy (we might expect a hard copy)" << std::endl;
    exit(1);
    M_vector = other;
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


};
