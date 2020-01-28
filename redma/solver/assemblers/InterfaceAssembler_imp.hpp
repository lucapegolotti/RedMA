namespace RedMA
{


template <class InVectorType, class InMatrixType>
Interface<InVectorType, InMatrixType>::
Interface(SHP(AssemblerType) assemblerFather, const unsigned int& indexFather,
          SHP(AssemblerType) assemblerChild, const unsigned int& indexChild,
          const unsigned int& interfaceID) :
  M_assemblerFather(assemblerFather),
  M_indexFather(indexFather),
  M_assemblerChild(assemblerChild),
  M_indexChild(indexChild),
  M_ID(interfaceID)
{

}

template <class InVectorType, class InMatrixType>
InterfaceAssembler<InVectorType, InMatrixType>::
InterfaceAssembler(const DataContainer& data,
                   const Interface<InVectorType, InMatrixType>& interface) :
  M_data(data),
  M_interface(interface)
{
    setup();
}

template <class InVectorType, class InMatrixType>
void
InterfaceAssembler<InVectorType, InMatrixType>::
setup()
{
    LifeV::LifeChrono chrono;
    chrono.start();

    printlog(YELLOW, "[InterfaceAssembler] initialize interface"
                     " assembler ...", M_data.getVerbose());

    buildCouplingMatrices();

    std::string msg = "done, in ";
    msg += std::to_string(chrono.diff());
    msg += " seconds\n";
    printlog(YELLOW, msg, M_data.getVerbose());
}

template <class InVectorType, class InMatrixType>
void
InterfaceAssembler<InVectorType, InMatrixType>::
addContributionRhs(BlockVector<BlockVector<InVectorType>>& rhs,
                   const BlockVector<BlockVector<InVectorType>>& sol,
                   const unsigned int& nPrimalBlocks)
{
    unsigned int fatherID = M_interface.M_indexFather;
    unsigned int childID = M_interface.M_indexChild;
    unsigned int interfaceID = M_interface.M_ID;
    SHP(aAssembler<InVectorType COMMA InMatrixType>) assemblerFather;
    SHP(aAssembler<InVectorType COMMA InMatrixType>) assemblerChild;
    assemblerFather = M_interface.M_assemblerFather;
    assemblerChild = M_interface.M_assemblerChild;

    // we have (-1) because we are solving H un+1 = F(.) and coupling is in F
    rhs.block(fatherID) -= M_fatherBT * sol.block(nPrimalBlocks + interfaceID);
    rhs.block(childID)  -= M_childBT * sol.block(nPrimalBlocks + interfaceID);
    assemblerFather->getBCManager()->apply0DirichletBCs(rhs.block(fatherID),
                                                        assemblerFather->getFESpaceBCs(),
                                                        assemblerFather->getComponentBCs());

    assemblerChild->getBCManager()->apply0DirichletBCs(rhs.block(childID),
                                                       assemblerChild->getFESpaceBCs(),
                                                       assemblerChild->getComponentBCs());

    rhs.block(nPrimalBlocks + interfaceID) -= M_fatherB * sol.block(fatherID);
    rhs.block(nPrimalBlocks + interfaceID) -= M_childB * sol.block(childID);
}

template <class InVectorType, class InMatrixType>
void
InterfaceAssembler<InVectorType, InMatrixType>::
addContributionJacobianRhs(BlockMatrix<BlockMatrix<InMatrixType>>& jac,
                           const BlockVector<BlockVector<InVectorType>>& sol,
                           const unsigned int& nPrimalBlocks)
{
    unsigned int fatherID = M_interface.M_indexFather;
    unsigned int childID = M_interface.M_indexChild;
    unsigned int interfaceID = M_interface.M_ID;

    // hard copy, otherwise we flip the sign of the matrices every time this
    // function is called
    jac.block(fatherID, nPrimalBlocks + interfaceID).hardCopy(M_fatherBT);
    jac.block(childID,  nPrimalBlocks + interfaceID).hardCopy(M_childBT);
    jac.block(nPrimalBlocks + interfaceID, fatherID).hardCopy(M_fatherB);
    jac.block(nPrimalBlocks + interfaceID,  childID).hardCopy(M_childBT);

    jac.block(fatherID, nPrimalBlocks + interfaceID) *= (-1);
    jac.block(childID,  nPrimalBlocks + interfaceID) *= (-1);
    jac.block(nPrimalBlocks + interfaceID, fatherID) *= (-1);
    jac.block(nPrimalBlocks + interfaceID,  childID) *= (-1);
}

}
