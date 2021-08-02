#include "DoubleVector.hpp"

namespace RedMA
{
DoubleVector::
DoubleVector()
{
    M_double = 0;
}

/*void
DoubleVector::
add(shp<aMatrix> other)
{
    if (!other->isZero())
    {
        double value = convert<DoubleVector>(other)->getValue();
        M_double += value;
    }
}*/

void
DoubleVector::
multiplyByScalar(const double& coeff)
{
    M_double *= coeff;
}

/*shp<aVector>
DoubleVector::
multiplyByVector(shp<aVector> vector)
{
    shp<DoubleVector> res(new DoubleVector());
    if (!vector->isZero())
    {
        double value = convert<DoubleVector>(vector)->getValue();

        res->setValue(M_double * value);
    }
    return res;
}*/

/*shp<aMatrix>
DoubleVector::
multiplyByMatrix(shp<aMatrix> other)
{
    shp<DoubleVector> res(new DoubleVector());
    if (!other->isZero())
    {
        double value = convert<DoubleVector>(other)->getValue();
        res->setValue(M_double * value);
    }
    return res;
}*/

void
DoubleVector::
dump(std::string namefile) const
{
    throw new Exception("Method dump not implemented for class DoubleVector!");
}

bool
DoubleVector::
isZero() const
{
    return false; //M_double == 0;
}

double
DoubleVector::
operator()(unsigned int index)
{
    return M_double;
}

void
DoubleVector::
add(shp<aVector> other)
{
    if (!other->isZero())
    {
        double value = convert<DoubleVector>(other)->getValue();
        M_double += value;
    }
}

void
DoubleVector::
shallowCopy(shp<aDataWrapper> other)
{
    double value = convert<DoubleVector>(other)->getValue();
    M_double = value;
}

void
DoubleVector::
deepCopy(shp<aDataWrapper> other)
{
    double value = 0;
    if (other)
        value = convert<DoubleVector>(other)->getValue();
    M_double = value;
}

std::string
DoubleVector::
getString(const char& delimiter) const
{
    throw new Exception("Method getString not implemented for class DoubleVector!");
}

DoubleVector*
DoubleVector::
clone() const
{
    DoubleVector* value = new DoubleVector();
    value->setValue(M_double);
    return value;
    // throw new Exception("Method clone not implemented for class DoubleVector!");
}

shp<void>
DoubleVector::
data() const
{
    shp<double> res(new double(M_double));
    return res;
}

// specification of template function to avoid ambiguous cast
template<>
shp<DoubleVector> convert(shp<aDataWrapper> container)
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
    return spcast<DoubleVector>(spcast<aVector>(container));
}

}
