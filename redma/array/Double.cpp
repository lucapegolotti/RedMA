#include "Double.hpp"

namespace RedMA
{
Double::
Double()
{
    M_double = 0;
}

void
Double::
add(std::shared_ptr<aMatrix> other)
{
    if (!other->isZero())
    {
        double value = convert<Double>(other)->getValue();
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
        double value = convert<Double>(vector)->getValue();

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
        double value = convert<Double>(other)->getValue();
        res->setValue(M_double * value);
    }
    return res;
}

void
Double::
dump(std::string namefile) const
{

}

bool
Double::
isZero() const
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
        double value = value = convert<Double>(other)->getValue();
        M_double += value;
    }
}

void
Double::
shallowCopy(std::shared_ptr<aDataWrapper> other)
{
    double value = convert<Double>(other)->getValue();
    M_double = value;
}

void
Double::
deepCopy(std::shared_ptr<aDataWrapper> other)
{
    double value = 0;
    if (other)
        value = convert<Double>(other)->getValue();
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

aDataWrapper*
Double::
clone() const
{
    throw new Exception("Method clone not implemented for Double!");
    // Double* retDouble = new Double();
    // retDouble->setValue(M_double);
    // return retDouble;
}

std::shared_ptr<void>
Double::
data() const
{
    std::shared_ptr<double> res(new double(M_double));
    return res;
}

// specification of template function to avoid ambiguous cast
template<>
shp<Double> convert(shp<aDataWrapper> container)
{
    if (container->type() != DOUBLE)
    {
        std::string msg = "Error in convert: converting ";
        msg += std::to_string(container->type());
        msg += " to ";
        msg += std::to_string(DOUBLE);
        msg += "\n";
        throw new Exception(msg);
    }
    return spcast<Double>(spcast<aVector>(container));
}

}
