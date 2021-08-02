#include "PressureDrop.hpp"

namespace RedMA
{

PressureDrop::
PressureDrop(const double& C, const double& Rd) :
  M_C(C),
  M_R(Rd)
{
}


shp<aVector>
PressureDrop::
getZeroVector() const
{
    shp<BlockVector> retVec(new BlockVector(1));
    shp<DoubleVector> value(new DoubleVector());
    value->setValue(0);
    retVec->setBlock(0, value);

    return retVec;
}

shp<aMatrix>
PressureDrop::
getMass(const double& time, const shp<aVector>& sol)
{
    shp<BlockMatrix> mass(new BlockMatrix(1,1));
    shp<DoubleMatrix> one(new DoubleMatrix());
    one->setValue(1);
    mass->setBlock(0,0,one);

    return mass;
}

shp<aMatrix>
PressureDrop::
getPressureMass(const double& time, const shp<aVector>& sol)
{
    return this->getMass(time, sol);
}

shp<aMatrix>
PressureDrop::
getMassJacobian(const double& time, const shp<aVector>& sol)
{
    shp<BlockMatrix> massJac(new BlockMatrix(1,1));
    return massJac;
}

shp<aVector>
PressureDrop::
getRightHandSide(const double& time, const shp<aVector>& sol)
{
    shp<BlockVector> retVec(new BlockVector(1));
    retVec->deepCopy(sol);
    retVec->block(0)->multiplyByScalar(-1.0 / (M_C * M_R));
    double v = spcast<DoubleVector>(retVec->block(0))->getValue();
    spcast<DoubleVector>(retVec->block(0))->setValue(v + M_Q / M_C);

    return retVec;
}

shp<aMatrix>
PressureDrop::
getJacobianRightHandSide(const double& time, const shp<aVector>& sol)
{
    shp<BlockMatrix> retMat(new BlockMatrix(1,1));
    shp<DoubleMatrix> value(new DoubleMatrix());
    value->setValue(-1.0 / (M_C * M_R));
    retMat->setBlock(0,0,value);

    return retMat;
}

}
