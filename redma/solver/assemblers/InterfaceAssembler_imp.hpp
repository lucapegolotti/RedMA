namespace RedMA
{

template <class InVectorType, class InMatrixType>
Interface<InVectorType, InMatrixType>::
Interface(SHP(AssemblerType) assemblerFather, const int& indexFather,
          SHP(AssemblerType) assemblerChild, const int& indexChild,
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

    M_stabilizationCoupling = M_data("coupling/stab_coefficient", 0.0);

    buildCouplingMatrices();

    std::string msg = "done, in ";
    msg += std::to_string(chrono.diff());
    msg += " seconds\n";
    printlog(YELLOW, msg, M_data.getVerbose());
}

template <class InVectorType, class InMatrixType>
void
InterfaceAssembler<InVectorType, InMatrixType>::
addContributionRhs(const double& time,
                   BlockVector<BlockVector<InVectorType>>& rhs,
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

    if (M_stabilizationCoupling > THRESHOLDSTAB)
    {
        rhs.block(nPrimalBlocks + interfaceID) -= (M_stabFather * sol.block(fatherID)) * (0.5 * M_stabilizationCoupling);
        rhs.block(nPrimalBlocks + interfaceID) -= (M_stabChild * sol.block(childID)) * (0.5 * M_stabilizationCoupling);

        rhs.block(nPrimalBlocks + interfaceID) -=
        sol.block(nPrimalBlocks + interfaceID) * M_stabilizationCoupling;
    }
}

template <class InVectorType, class InMatrixType>
double
InterfaceAssembler<InVectorType, InMatrixType>::
checkStabilizationTerm(const BlockVector<BlockVector<InVectorType>>& sol,
                       const unsigned int& nPrimalBlocks)
{
    if (M_stabilizationCoupling > THRESHOLDSTAB &&
        M_interface.M_assemblerFather && M_interface.M_assemblerChild)
    {
        unsigned int fatherID = M_interface.M_indexFather;
        unsigned int childID = M_interface.M_indexChild;
        unsigned int interfaceID = M_interface.M_ID;

        BlockVector<BlockVector<InVectorType>> res;
        res.resize(1);

        res.block(0) -= (M_stabFather * sol.block(fatherID)) * 0.5;
        res.block(0) -= (M_stabChild * sol.block(childID)) * 0.5;

        // std::cout << "--------" << std::endl << std::flush;
        // res.block(0).block(0).data()->showMe();
        // std::cout << "++++++++" << std::endl << std::flush;
        // sol.block(nPrimalBlocks + interfaceID).block(0).data()->showMe();
        //
        // std::cout << "stab term stress " << res.norm2() << std::endl << std::flush;
        // std::cout << "stab term lagrange " << sol.block(nPrimalBlocks + interfaceID).norm2() << std::endl << std::flush;

        res.block(0) -= sol.block(nPrimalBlocks + interfaceID);

        std::string msg = "[InterfaceAssembler] interface ID = ";
        msg += std::to_string(interfaceID);
        msg += ", stab term norm = ";
        msg += std::to_string(res.norm2());
        msg += "\n";
        printlog(MAGENTA, msg, M_data.getVerbose());
    }
}

template <class InVectorType, class InMatrixType>
void
InterfaceAssembler<InVectorType, InMatrixType>::
addContributionJacobianRhs(const double& time,
                           BlockMatrix<BlockMatrix<InMatrixType>>& jac,
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
    jac.block(nPrimalBlocks + interfaceID,  childID).hardCopy(M_childB);

    jac.block(fatherID, nPrimalBlocks + interfaceID) *= (-1);
    jac.block(childID,  nPrimalBlocks + interfaceID) *= (-1);
    jac.block(nPrimalBlocks + interfaceID, fatherID) *= (-1);
    jac.block(nPrimalBlocks + interfaceID,  childID) *= (-1);

    if (M_stabilizationCoupling > THRESHOLDSTAB)
    {
        jac.block(nPrimalBlocks + interfaceID, fatherID) += (M_stabFather * (-0.5 * M_stabilizationCoupling));
        jac.block(nPrimalBlocks + interfaceID,  childID) += (M_stabChild * (-0.5 * M_stabilizationCoupling));

        jac.block(nPrimalBlocks + interfaceID, nPrimalBlocks + interfaceID).hardCopy(M_identity * (-1.0 * M_stabilizationCoupling));
    }
}

}
