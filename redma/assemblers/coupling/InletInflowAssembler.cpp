#include "InletInflowAssembler.hpp"

namespace RedMA
{

InletInflowAssembler::
InletInflowAssembler(const DataContainer& data,
                     const Interface& interface) :
  InterfaceAssembler(data, interface)
{
}

void
InletInflowAssembler::
addContributionJacobianRhs(const double& time,
                           shp<BlockMatrix> jac,
                           shp<BlockVector> sol,
                           const unsigned int& nPrimalBlocks)
{
    unsigned int fatherID = this->M_interface.M_indexFather;
    unsigned int childID = this->M_interface.M_indexChild;
    unsigned int interfaceID = this->M_interface.M_ID;

    // hard copy, otherwise we flip the sign of the matrices every time this
    // function is called
    jac->setBlock(childID,  nPrimalBlocks + interfaceID,shp<BlockMatrix>(this->M_childBT->clone()));
    jac->setBlock(nPrimalBlocks + interfaceID,  childID,shp<BlockMatrix>(this->M_childB->clone()));

    jac->block(childID,  nPrimalBlocks + interfaceID)->multiplyByScalar(-1);
    jac->block(nPrimalBlocks + interfaceID,  childID)->multiplyByScalar(-1);

    // if (this->M_stabilizationCoupling > THRESHOLDSTAB)
    // {
    //     jac.block(nPrimalBlocks + interfaceID,  childID) += (this->M_stabChild * (-1.0 * this->M_stabilizationCoupling));
    //
    //     jac.block(nPrimalBlocks + interfaceID, nPrimalBlocks + interfaceID).deepCopy(this->M_identity * (-1.0 * this->M_stabilizationCoupling));
    // }
}

void
InletInflowAssembler::
addContributionRhs(const double& time, shp<BlockVector> rhs, shp<BlockVector> sol,
                   const unsigned int& nPrimalBlocks)
{
    unsigned int childID = this->M_interface.M_indexChild;
    unsigned int interfaceID = this->M_interface.M_ID;
    shp<aAssembler> assemblerChild = this->M_interface.M_assemblerChild;

    auto temp = this->M_childBT->multiplyByVector(sol->block(nPrimalBlocks + interfaceID));
    temp->multiplyByScalar(-1);
    rhs->block(childID)->add(temp);

    temp = this->M_childB->multiplyByVector(sol->block(childID));
    temp->multiplyByScalar(-1);
    if (rhs->block(nPrimalBlocks + interfaceID)->isZero())
        rhs->setBlock(nPrimalBlocks + interfaceID,temp);
    else
        rhs->block(nPrimalBlocks + interfaceID)->add(temp);

    // rhs->block(nPrimalBlocks + interfaceID) -= this->M_childB * sol.block(childID);
    temp = this->M_childBfe->multiplyByVector(assemblerChild->getLifting(time));
    temp->multiplyByScalar(1.0/M_data.getInflow()(time));
    temp->dump("couplingRHS");

    if (assemblerChild->getRBBases())
        temp = assemblerChild->getRBBases()->projectOnLagrangeSpace(spcast<BlockVector>(temp));
    rhs->block(nPrimalBlocks + interfaceID)->add(temp);

    // // here + because the stabilization term is (stress - lagrange) => hence
    // // + lagrange at rhs
    // if (this->M_stabilizationCoupling > THRESHOLDSTAB)
    // {
    //     rhs.block(nPrimalBlocks + interfaceID) -= (this->M_stabChild * sol.block(childID)) * (1.0 * this->M_stabilizationCoupling);
    //
    //     rhs.block(nPrimalBlocks + interfaceID) -=
    //     sol.block(nPrimalBlocks + interfaceID) * this->M_stabilizationCoupling;
    // }
}

}
