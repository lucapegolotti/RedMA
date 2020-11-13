#include "PressureDrop.hpp"

namespace RedMA
{

PressureDrop::
PressureDrop(const double& C, const double& Rp, const double& Rd) :
  M_C(C),
  M_Rp(Rp),
  M_Rd(Rd)
{

}


SHP(aVector)
PressureDrop::
getZeroVector() const
{
    SHP(BlockVector) retVec(new BlockVector(1));
    SHP(Double) value(new Double());
    value->setValue(0);
    retVec->setBlock(0,value);
    return retVec;
}

SHP(aMatrix)
PressureDrop::
getMass(const double& time, const SHP(aVector)& sol)
{
    SHP(BlockMatrix) mass(new BlockMatrix(1,1));
    SHP(Double) one(new Double());
    one->setValue(1.0);
    mass->setBlock(0,0,one);

    return mass;
}

SHP(aMatrix)
PressureDrop::
getMassJacobian(const double& time, const SHP(aVector)& sol)
{
    SHP(BlockMatrix) massJac(new BlockMatrix(1,1));
    return massJac;
}

SHP(aVector)
PressureDrop::
getRightHandSide(const double& time, const SHP(aVector)& sol)
{
    SHP(BlockVector) retVec(new BlockVector(1));
    retVec->hardCopy(sol);
    std::static_pointer_cast<Double>(retVec->block(0))->multiplyByScalar(-1.0 / (M_C * M_Rd));
    double v = std::static_pointer_cast<Double>(retVec->block(0))->getValue();
    // std::static_pointer_cast<Double>(retVec->block(0))->setValue(v + M_Q / M_C);
    return retVec;
}

SHP(aMatrix)
PressureDrop::
getJacobianRightHandSide(const double& time, const SHP(aVector)& sol)
{
    SHP(BlockMatrix) retMat(new BlockMatrix(1,1));
    SHP(Double) value(new Double());
    // value->setValue(-1.0 / (M_C * M_Rd));
    retMat->setBlock(0,0,value);
    return retMat;
}

}
