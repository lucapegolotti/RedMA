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


BlockVector
PressureDrop::
getZeroVector() const
{
    // BlockVector<Double> retVec(1);
    // retVec.block(0).data() = 0;
    // return retVec;
}

BlockMatrix
PressureDrop::
getMass(const double& time, const BlockVector& sol)
{
    // BlockMatrix<Double> mass(1,1);
    // mass.block(0,0) = 1.0;
    // return mass;
}

BlockMatrix
PressureDrop::
getMassJacobian(const double& time, const BlockVector& sol)
{
    // BlockMatrix<Double> massJac(1,1);
    // return massJac;
}

BlockVector
PressureDrop::
getRightHandSide(const double& time, const BlockVector& sol)
{
    // BlockVector retVec;
    // retVec.hardCopy(sol);
    // retVec.block(0).data() *= (-1.0 / (M_C * M_Rd));
    // retVec.block(0).data() += M_Q / M_C;
    // return retVec;
}

BlockMatrix
PressureDrop::
getJacobianRightHandSide(const double& time, const BlockVector& sol)
{
    // BlockMatrix retMat(1,1);
    // retMat.block(0,0).data() = -1.0 / (M_C * M_Rd);
    // return retMat;
}

}