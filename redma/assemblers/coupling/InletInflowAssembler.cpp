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
                           BlockMatrix& jac,
                           const BlockVector& sol,
                           const unsigned int& nPrimalBlocks)
{
    // unsigned int fatherID = this->M_interface.M_indexFather;
    // unsigned int childID = this->M_interface.M_indexChild;
    // unsigned int interfaceID = this->M_interface.M_ID;
    //
    // // hard copy, otherwise we flip the sign of the matrices every time this
    // // function is called
    // jac.block(childID,  nPrimalBlocks + interfaceID).hardCopy(this->M_childBT);
    // jac.block(nPrimalBlocks + interfaceID,  childID).hardCopy(this->M_childB);
    //
    // jac.block(childID,  nPrimalBlocks + interfaceID) *= (-1);
    // jac.block(nPrimalBlocks + interfaceID,  childID) *= (-1);
    //
    // if (this->M_stabilizationCoupling > THRESHOLDSTAB)
    // {
    //     jac.block(nPrimalBlocks + interfaceID,  childID) += (this->M_stabChild * (-1.0 * this->M_stabilizationCoupling));
    //
    //     jac.block(nPrimalBlocks + interfaceID, nPrimalBlocks + interfaceID).hardCopy(this->M_identity * (-1.0 * this->M_stabilizationCoupling));
    // }
}

}
