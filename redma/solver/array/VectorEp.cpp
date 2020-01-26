#include "VectorEp.hpp"

namespace RedMA
{

VectorEp::
VectorEp()
{

}

VectorEp
VectorEp::
operator+(const VectorEp& other)
{

}

VectorEp&
VectorEp::
operator+=(const VectorEp& other)
{

}

VectorEp&
VectorEp::
operator-=(const VectorEp& other)
{

}

VectorEp&
VectorEp::
operator=(const std::shared_ptr<VECTOREPETRA>& other)
{
    M_vector = other;
}

VectorEp&
VectorEp::
operator*=(const double& coeff)
{

}

void
VectorEp::
hardCopy(const VectorEp& other)
{

}

void
VectorEp::
softCopy(const VectorEp& other)
{
    
}

double
VectorEp::
norm2() const
{

}

std::shared_ptr<VECTOREPETRA>&
VectorEp::
data()
{
    return M_vector;
}

};
