#include "WindkesselPressureDrop.hpp"

namespace RedMA
{

WindkesselPressureDrop::
WindkesselPressureDrop(const double& C, const double& Rd) :
  M_C(C),
  M_R(Rd)
{
}

shp<aVector>
WindkesselPressureDrop::
getZeroVector() const
{
    shp<BlockVector> retVec(new BlockVector(1));
    shp<DoubleVector> value(new DoubleVector());
    value->setValue(0);
    retVec->setBlock(0, value);

    return retVec;
}

shp<aMatrix>
WindkesselPressureDrop::
getMass(const double& time, const shp<aVector>& sol)
{
    shp<BlockMatrix> mass(new BlockMatrix(1,1));
    shp<DoubleMatrix> value(new DoubleMatrix());
    value->setValue(1.0);
    mass->setBlock(0,0,value);

    return mass;
}

shp<aMatrix>
WindkesselPressureDrop::
getMassJacobian(const double& time, const shp<aVector>& sol)
{
    shp<BlockMatrix> massJac(new BlockMatrix(1,1));
    return massJac;
}

shp<aVector>
WindkesselPressureDrop::
getRightHandSide(const double& time, const shp<aVector>& sol)
{
    shp<BlockVector> retVec(new BlockVector(1));
    retVec->deepCopy(sol);
    retVec->block(0)->multiplyByScalar(-1.0 / (M_R * M_C));
    double v = spcast<DoubleVector>(retVec->block(0))->getValue();
    spcast<DoubleVector>(retVec->block(0))->setValue(v + M_Q / M_C);

    return retVec;
}

shp<aMatrix>
WindkesselPressureDrop::
getJacobianRightHandSide(const double& time, const shp<aVector>& sol)
{
    shp<BlockMatrix> retMat(new BlockMatrix(1,1));
    shp<DoubleMatrix> value(new DoubleMatrix());
    value->setValue(-1.0 / (M_R * M_C));
    retMat->setBlock(0,0,value);

    return retMat;
}

}
