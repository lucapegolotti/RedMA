namespace RedMA
{

template <class InVectorType, class InMatrixType>
InletInflowAssembler<InVectorType, InMatrixType>::
InletInflowAssembler(const DataContainer& data,
                     const Interface<InVectorType, InMatrixType>& interface) :
  InterfaceAssembler<InVectorType, InMatrixType>(data, interface)
{
}

template <class InVectorType, class InMatrixType>
void
InletInflowAssembler<InVectorType, InMatrixType>::
addContributionRhs(const double& time,
                   BlockVector<BlockVector<InVectorType>>& rhs,
                   const BlockVector<BlockVector<InVectorType>>& sol,
                   const unsigned int& nPrimalBlocks)
{
    unsigned int childID = this->M_interface.M_indexChild;
    unsigned int interfaceID = this->M_interface.M_ID;
    SHP(aAssembler<InVectorType COMMA InMatrixType>) assemblerChild;
    assemblerChild = this->M_interface.M_assemblerChild;

    rhs.block(childID)  -= this->M_childBT * sol.block(nPrimalBlocks + interfaceID);

    assemblerChild->getBCManager()->apply0DirichletBCs(rhs.block(childID),
                                                       assemblerChild->getFESpaceBCs(),
                                                       assemblerChild->getComponentBCs());

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

template <class InVectorType, class InMatrixType>
void
InletInflowAssembler<InVectorType, InMatrixType>::
addContributionJacobianRhs(const double& time,
                           BlockMatrix<BlockMatrix<InMatrixType>>& jac,
                           const BlockVector<BlockVector<InVectorType>>& sol,
                           const unsigned int& nPrimalBlocks)
{
    unsigned int fatherID = this->M_interface.M_indexFather;
    unsigned int childID = this->M_interface.M_indexChild;
    unsigned int interfaceID = this->M_interface.M_ID;

    // hard copy, otherwise we flip the sign of the matrices every time this
    // function is called
    jac.block(childID,  nPrimalBlocks + interfaceID).hardCopy(this->M_childBT);
    jac.block(nPrimalBlocks + interfaceID,  childID).hardCopy(this->M_childB);

    jac.block(childID,  nPrimalBlocks + interfaceID) *= (-1);
    jac.block(nPrimalBlocks + interfaceID,  childID) *= (-1);

    if (this->M_stabilizationCoupling > THRESHOLDSTAB)
    {
        jac.block(nPrimalBlocks + interfaceID,  childID) += (this->M_stabChild * (-1.0 * this->M_stabilizationCoupling));

        jac.block(nPrimalBlocks + interfaceID, nPrimalBlocks + interfaceID).hardCopy(this->M_identity * (-1.0 * this->M_stabilizationCoupling));
    }
}

}