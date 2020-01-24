#include "MatrixEp.hpp"

namespace RedMA
{

MatrixEp::
MatrixEp()
{

}

MatrixEp
MatrixEp::
operator+(const MatrixEp& other)
{

}

MatrixEp&
MatrixEp::
operator+=(const MatrixEp& other)
{

}

MatrixEp&
MatrixEp::
operator*=(const double& coeff)
{

}

void
MatrixEp::
hardCopy(const MatrixEp& other)
{

}

void
MatrixEp::
softCopy(const MatrixEp& other)
{

}

VectorEp
MatrixEp::
operator*(const VectorEp& vector)
{

}

void
MatrixEp::
getRowProperty(std::shared_ptr<LifeV::MapEpetra>& outMap)
{

}

void
MatrixEp::
getColProperty(std::shared_ptr<LifeV::MapEpetra>& outMap)
{

}

std::shared_ptr<MATRIXEPETRA>&
MatrixEp::
data()
{
    return M_matrix;
}

}
