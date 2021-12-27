#include "CoronaryPressureDrop.hpp"

namespace RedMA
{

CoronaryPressureDrop::
CoronaryPressureDrop(const double &Ca, const double &Cim,
                     const double &Ram, const double &Rvm,
                     const double &Rv) :
    M_Ca(Ca), M_Cim(Cim), M_Ram(Ram), M_Rvm(Rvm), M_Rv(Rv)
{
}

shp<aVector>
CoronaryPressureDrop::
getZeroVector() const
{
    shp<BlockVector> retVec(new BlockVector(2));

    shp<DoubleVector> value0(new DoubleVector());
    value0->setValue(0);
    retVec->setBlock(0, value0);

    shp<DoubleVector> value1(new DoubleVector());
    value1->setValue(0);
    retVec->setBlock(0, value1);

    return retVec;
}


shp<aMatrix>
CoronaryPressureDrop::
getMass(const double &time, const shp<aVector>& sol)
{
    shp<BlockMatrix> mass(new BlockMatrix(2,2));

    shp<DoubleMatrix> value0(new DoubleMatrix());
    value0->setValue(M_Ca);
    mass->setBlock(0,0,value0);

    shp<DoubleMatrix> value1(new DoubleMatrix());
    value1->setValue(M_Cim);
    mass->setBlock(1,1,value1);

    return mass;
}

shp<aMatrix>
CoronaryPressureDrop::
getMassJacobian(const double& time, const shp<aVector>& sol)
{
    shp<BlockMatrix> massJac(new BlockMatrix(2,2));
    return massJac;
}

shp<aVector>
CoronaryPressureDrop::
getRightHandSide(const double& time, const shp<aVector>& sol)
{
    shp<BlockVector> retVec(new BlockVector(2));

    double Pi1 = spcast<DoubleVector>(convert<BlockVector>(sol)->block(0))->getValue();
    double Pi2 = spcast<DoubleVector>(convert<BlockVector>(sol)->block(1))->getValue();

    shp<DoubleVector> value0(new DoubleVector());
    value0->setValue(-Pi1 / M_Ram +
                     Pi2 / M_Ram + M_Q +
                     M_Pim / M_Ram);
    retVec->setBlock(0, value0);

    shp<DoubleVector> value1(new DoubleVector());
    value1->setValue(Pi1 / M_Ram -
                     Pi2 * (1.0 / M_Ram + 1.0 / (M_Rvm + M_Rv)) -
                     M_Pim * (1.0 / M_Ram + 1.0 / (M_Rvm + M_Rv)));
    retVec->setBlock(1, value1);
    
    return retVec;
}

shp<aMatrix>
CoronaryPressureDrop::
getJacobianRightHandSide(const double& time, const shp<aVector>& sol)
{
    shp<BlockMatrix> retMat(new BlockMatrix(2,2));
    
    shp<DoubleMatrix> value00(new DoubleMatrix());
    value00->setValue(-1.0 / M_Ram);
    retMat->setBlock(0,0,value00);

    shp<DoubleMatrix> value01(new DoubleMatrix());
    value01->setValue(1.0 / M_Ram);
    retMat->setBlock(0,1,value01);

    shp<DoubleMatrix> value10(new DoubleMatrix());
    value10->setValue(1.0 / M_Ram);
    retMat->setBlock(1,0,value10);

    shp<DoubleMatrix> value11(new DoubleMatrix());
    value11->setValue(- 1.0 / M_Ram - 1.0 / (M_Rvm + M_Rv));
    retMat->setBlock(1,1,value11);

    return retMat;
}


}  // namespace RedMA
