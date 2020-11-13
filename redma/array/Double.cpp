#include "Double.hpp"

namespace RedMA
{
Double::
Double() :
  aMatrix(DOUBLE),
  aVector(DOUBLE)
{
    M_double = 0;
}

void
Double::
add(std::shared_ptr<aMatrix> other)
{
    if (!other->isZero())
    {
        aMatrix::checkType(other, DOUBLE);
        double value = std::static_pointer_cast<Double>(other)->getValue();
        M_double += value;
    }
}

void
Double::
multiplyByScalar(const double& coeff)
{
    M_double *= coeff;
}

std::shared_ptr<aVector>
Double::
multiplyByVector(std::shared_ptr<aVector> vector)
{
    std::shared_ptr<Double> res(new Double());
    if (!vector->isZero())
    {
        aMatrix::checkType(vector, DOUBLE);
        double value = std::static_pointer_cast<Double>(vector)->getValue();

        res->setValue(M_double * value);
    }
    return res;
}

std::shared_ptr<aMatrix>
Double::
multiplyByMatrix(std::shared_ptr<aMatrix> other)
{
    std::shared_ptr<Double> res(new Double());
    if (!other->isZero())
    {
        aMatrix::checkType(other, DOUBLE);
        double value = std::static_pointer_cast<Double>(other)->getValue();
        res->setValue(M_double * value);
    }
    return res;
}

void
Double::
dump(std::string namefile) const
{

}

void
Double::
softCopy(std::shared_ptr<aMatrix> other)
{
    aMatrix::checkType(other, DOUBLE);
    double value = std::static_pointer_cast<Double>(other)->getValue();
    M_double = value;
}

void
Double::
hardCopy(std::shared_ptr<aMatrix> other)
{
    std::cout << "hardCopy" << std::endl << std::flush;
    aMatrix::checkType(other, DOUBLE);
    double value = 0;
    if (other)
        value = std::static_pointer_cast<Double>(other)->getValue();
    M_double = value;
}

bool
Double::
isZero()
{
    return M_double == 0;
}

double
Double::
operator()(unsigned int index)
{
    return M_double;
}

void
Double::
add(std::shared_ptr<aVector> other)
{
    if (!other->isZero())
    {
        aVector::checkType(other, DOUBLE);
        double value = std::static_pointer_cast<Double>(other)->getValue();
        M_double += value;
    }
}

void
Double::
softCopy(std::shared_ptr<aVector> other)
{
    aVector::checkType(other, DOUBLE);
    double value = std::static_pointer_cast<Double>(other)->getValue();
    M_double = value;
}

void
Double::
hardCopy(std::shared_ptr<aVector> other)
{
    aVector::checkType(other, DOUBLE);
    double value = 0;
    if (other)
        value = std::static_pointer_cast<Double>(other)->getValue();
    M_double = value;
}

std::string
Double::
getString(const char& delimiter) const
{

}

double
Double::
norm2() const
{
    return M_double > 0 ? M_double : -M_double;
}

aVector*
Double::
cloneVector() const
{
    Double* retVector = new Double();
    retVector->setValue(M_double);
    return retVector;
}

aMatrix*
Double::
clone() const
{
    Double* retMatrix = new Double();
    retMatrix->setValue(M_double);
    return retMatrix;
}

std::shared_ptr<void>
Double::
data() const
{
    std::shared_ptr<double> res(new double(M_double));
    return res;
}

}
