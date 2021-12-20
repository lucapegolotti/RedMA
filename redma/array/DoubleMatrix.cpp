#include "DoubleMatrix.hpp"

namespace RedMA
{
DoubleMatrix::
DoubleMatrix()
{
    M_double = 0;
}

void
DoubleMatrix::
add(shp<aMatrix> other)
{
    if (!other->isZero())
    {
        double value = convert<DoubleMatrix>(other)->getValue();
        M_double += value;
    }
}

void
DoubleMatrix::
multiplyByScalar(const double& coeff)
{
    M_double *= coeff;
}

shp<aVector>
DoubleMatrix::
multiplyByVector(shp<aVector> vector)
{
    shp<DoubleVector> res(new DoubleVector());
    if (!vector->isZero())
    {
        double value = convert<DoubleVector>(vector)->getValue();

        res->setValue(M_double * value);
    }
    return res;
}

shp<aMatrix>
DoubleMatrix::
multiplyByMatrix(shp<aMatrix> other)
{
    shp<DoubleMatrix> res(new DoubleMatrix());
    if (!other->isZero())
    {
        double value = convert<DoubleMatrix>(other)->getValue();
        res->setValue(M_double * value);
    }
    return res;
}

void
DoubleMatrix::
dump(std::string namefile) const
{
    throw new Exception("Method dump not implemented for class DoubleMatrix!");
}

bool
DoubleMatrix::
isZero() const
{
    return false; // M_double == 0;
}

void
DoubleMatrix::
shallowCopy(shp<aDataWrapper> other)
{
    double value = convert<DoubleMatrix>(other)->getValue();
    M_double = value;
}

void
DoubleMatrix::
deepCopy(shp<aDataWrapper> other)
{
    double value = 0;
    if (other)
        value = convert<DoubleMatrix>(other)->getValue();
    M_double = value;
}

DoubleMatrix*
DoubleMatrix::
clone() const
{
    DoubleMatrix* value = new DoubleMatrix();
    value->setValue(M_double);
    return value;
    // throw new Exception("Method clone not implemented for class DoubleMatrix!");
}

shp<void>
DoubleMatrix::
data() const
{
    shp<double> res(new double(M_double));
    return res;
}

    
} // namespace RedMA