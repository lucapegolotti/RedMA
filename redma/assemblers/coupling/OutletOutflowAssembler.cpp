#include "OutletOutflowAssembler.hpp"

namespace RedMA
{

OutletOutflowAssembler::
OutletOutflowAssembler(const DataContainer& data,
                     const Interface& interface,
                     const bool& addNoSlipBC) :
 InterfaceAssembler(data, interface, addNoSlipBC)
 {
 }

 void
 OutletOutflowAssembler::
 addContributionJacobianRhs(const double& time,
                            shp<BlockMatrix> jac,
                            shp<BlockVector> sol,
                            const unsigned int& nPrimalBlocks)
{
    unsigned int fatherID = this->M_interface.M_indexFather;
    unsigned int interfaceID = this->M_interface.M_ID;

    // hard copy, otherwise we flip the sign of the matrices every time this
    // function is called
    jac->setBlock(fatherID,  nPrimalBlocks + interfaceID,shp<BlockMatrix>(this->M_fatherBT->clone()));
    jac->setBlock(nPrimalBlocks + interfaceID,  fatherID,shp<BlockMatrix>(this->M_fatherB->clone()));

    jac->block(fatherID,  nPrimalBlocks + interfaceID)->multiplyByScalar(-1);
    jac->block(nPrimalBlocks + interfaceID,  fatherID)->multiplyByScalar(-1);
}

void
OutletOutflowAssembler::
addContributionRhs(const double& time, shp<BlockVector> rhs, shp<BlockVector> sol,
                   const unsigned int& nPrimalBlocks)
{
    unsigned int fatherID = this->M_interface.M_indexFather;
    unsigned int interfaceID = this->M_interface.M_ID;
    shp<aAssembler> assemblerFather = this->M_interface.M_assemblerFather;

    auto temp = this->M_fatherBT->multiplyByVector(sol->block(nPrimalBlocks + interfaceID));
    temp->multiplyByScalar(-1);
    rhs->block(fatherID)->add(temp);

    temp = this->M_fatherB->multiplyByVector(sol->block(fatherID));
    temp->multiplyByScalar(-1);
    if (rhs->block(nPrimalBlocks + interfaceID)->isZero())
        rhs->setBlock(nPrimalBlocks + interfaceID,temp);
    else
        rhs->block(nPrimalBlocks + interfaceID)->add(temp);
    
    temp = this->M_fatherBfe->multiplyByVector(assemblerFather->getLifting(time));
    temp->multiplyByScalar(-1); // correcting the sign
    if (assemblerFather->getRBBases())
        temp = assemblerFather->getRBBases()->projectOnLagrangeSpace(spcast<BlockVector>(temp));
    rhs->block(nPrimalBlocks + interfaceID)->add(temp);
}

}

