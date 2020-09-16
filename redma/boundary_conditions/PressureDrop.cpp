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
    // BlockVector<Double> retVec(1);
    // retVec.block(0).data() = 0;
    // return retVec;
}

SHP(aMatrix)
PressureDrop::
getMass(const double& time, const SHP(aVector)& sol)
{
    // BlockMatrix<Double> mass(1,1);
    // mass.block(0,0) = 1.0;
    // return mass;
}

SHP(aMatrix)
PressureDrop::
getMassJacobian(const double& time, const SHP(aVector)& sol)
{
    // BlockMatrix<Double> massJac(1,1);
    // return massJac;
}

SHP(aVector)
PressureDrop::
getRightHandSide(const double& time, const SHP(aVector)& sol)
{
    // SHP(BlockVector) retVec;
    // retVec.hardCopy(sol);
    // retVec.block(0).data() *= (-1.0 / (M_C * M_Rd));
    // retVec.block(0).data() += M_Q / M_C;
    // return retVec;
}

SHP(aMatrix)
PressureDrop::
getJacobianRightHandSide(const double& time, const SHP(aVector)& sol)
{
    // BlockMatrix retMat(1,1);
    // retMat.block(0,0).data() = -1.0 / (M_C * M_Rd);
    // return retMat;
}

}
