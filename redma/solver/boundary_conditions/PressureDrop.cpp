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


BlockVector<Double>
PressureDrop::
getZeroVector() const
{
    BlockVector<Double> retVec(1);
    retVec.block(0).data() = 0;
    return retVec;
}

BlockMatrix<Double>
PressureDrop::
getMass(const double& time, const BlockVector<Double>& sol)
{
    BlockMatrix<Double> mass(1,1);
    mass.block(0,0) = 1.0;
    return mass;
}

BlockMatrix<Double>
PressureDrop::
getMassJacobian(const double& time, const BlockVector<Double>& sol)
{
    BlockMatrix<Double> massJac(1,1);
    return massJac;
}

BlockVector<Double>
PressureDrop::
getRightHandSide(const double& time, const BlockVector<Double>& sol)
{
    BlockVector<Double> retVec;
    retVec.hardCopy(sol);
    retVec.block(0).data() *= (-1.0 / (M_C * M_Rd));
    retVec.block(0).data() += M_Q / M_C;
    return retVec;
}

BlockMatrix<Double>
PressureDrop::
getJacobianRightHandSide(const double& time, const BlockVector<Double>& sol)
{
    BlockMatrix<Double> retMat(1,1);
    retMat.block(0,0).data() = -1.0 / (M_C * M_Rd);
    return retMat;
}

}