#include "InletInflowAssembler.hpp"

namespace RedMA
{

template <>
void
InletInflowAssembler<VectorEp, MatrixEp>::
addContributionRhs(const double& time,
                   BlockVector<BlockVector<VectorEp>>& rhs,
                   const BlockVector<BlockVector<VectorEp>>& sol,
                   const unsigned int& nPrimalBlocks)
{
    unsigned int childID = this->M_interface.M_indexChild;
    unsigned int interfaceID = this->M_interface.M_ID;
    SHP(aAssembler<VectorEp COMMA MatrixEp>) assemblerChild;
    assemblerChild = this->M_interface.M_assemblerChild;

    rhs.block(childID)  -= this->M_childBT * sol.block(nPrimalBlocks + interfaceID);

    rhs.block(nPrimalBlocks + interfaceID) -= this->M_childB * sol.block(childID);
    rhs.block(nPrimalBlocks + interfaceID) += this->M_childB *
                                              assemblerChild->getLifting(time);
    // here + because the stabilization term is (stress - lagrange) => hence
    // + lagrange at rhs
    if (this->M_stabilizationCoupling > THRESHOLDSTAB)
    {
        rhs.block(nPrimalBlocks + interfaceID) -= (this->M_stabChild * sol.block(childID)) * (1.0 * this->M_stabilizationCoupling);

        rhs.block(nPrimalBlocks + interfaceID) -=
        sol.block(nPrimalBlocks + interfaceID) * this->M_stabilizationCoupling;
    }
}

}
