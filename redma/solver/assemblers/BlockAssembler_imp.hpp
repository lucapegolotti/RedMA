namespace RedMA
{

template <class InVectorType, class InMatrixType>
BlockAssembler<InVectorType, InMatrixType>::
BlockAssembler(const DataContainer& data, const TreeStructure& tree) :
  aAssembler<InVectorType, InMatrixType>(data),
  M_tree(tree)
{
    setup();
}

template <class InVectorType, class InMatrixType>
void
BlockAssembler<InVectorType, InMatrixType>::
checkStabTerm(const BlockVector<InVectorType>& sol) const
{
    for (auto as: M_dualAssemblers)
        as->checkStabilizationTerm(sol, M_primalAssemblers.size());
}

template <class InVectorType, class InMatrixType>
BlockVector<InVectorType>
BlockAssembler<InVectorType, InMatrixType>::
getLifting(const double& time) const
{
    BlockVector<InVectorType> retVec(M_numberBlocks);

    for (auto as : M_primalAssemblers)
        retVec.block(as.first) = as.second->getLifting(time);

    return retVec;
}

template <class InVectorType, class InMatrixType>
void
BlockAssembler<InVectorType, InMatrixType>::
apply0DirichletBCsMatrix(BlockMatrix<InMatrixType>& matrix, double diagCoeff) const
{
    for (auto as : M_primalAssemblers)
        as.second->apply0DirichletBCsMatrix(matrix.block(as.first, as.first), diagCoeff);
}

template <class InVectorType, class InMatrixType>
void
BlockAssembler<InVectorType, InMatrixType>::
apply0DirichletBCs(BlockVector<InVectorType>& initialGuess) const
{
    for (auto as : M_primalAssemblers)
        as.second->apply0DirichletBCs(initialGuess.block(as.first));
}

template <class InVectorType, class InMatrixType>
void
BlockAssembler<InVectorType, InMatrixType>::
applyDirichletBCs(const double& time, BlockVector<InVectorType>& initialGuess) const
{
    for (auto as : M_primalAssemblers)
        as.second->applyDirichletBCs(time, initialGuess.block(as.first));
}

template <class InVectorType, class InMatrixType>
BlockVector<InVectorType>
BlockAssembler<InVectorType, InMatrixType>::
getZeroVector() const
{
    BlockVector<InVectorType> retVec;
    retVec.resize(M_numberBlocks);

    for (auto as : M_primalAssemblers)
        retVec.block(as.first).softCopy(as.second->getZeroVector());

    unsigned int count = M_primalAssemblers.size();
    for (auto as : M_dualAssemblers)
    {
        retVec.block(count).softCopy(as->getZeroVector());
        count++;
    }

    return retVec;
}


template <class InVectorType, class InMatrixType>
void
BlockAssembler<InVectorType, InMatrixType>::
exportSolution(const double& t, const BlockVector<InVectorType>& sol)
{
    for (auto as : M_primalAssemblers)
        as.second->exportSolution(t, sol.block(as.first));
}

template <class InVectorType, class InMatrixType>
void
BlockAssembler<InVectorType, InMatrixType>::
setExtrapolatedSolution(const BlockVector<InVectorType>& exSol)
{
    for (auto as : M_primalAssemblers)
        as.second->setExtrapolatedSolution(exSol.block(as.first));
}

template <class InVectorType, class InMatrixType>
void
BlockAssembler<InVectorType, InMatrixType>::
postProcess(const double& t, const BlockVector<InVectorType>& sol)
{
    for (auto as : M_primalAssemblers)
        as.second->postProcess(t, sol.block(as.first));

    if (this->M_data("coupling/check_stabterm", false))
        checkStabTerm(sol);

}

template <class InVectorType, class InMatrixType>
BlockMatrix<InMatrixType>
BlockAssembler<InVectorType, InMatrixType>::
getMass(const double& time, const BlockVector<InVectorType>& sol)
{
    BlockMatrix<InMatrixType> mass;
    mass.resize(M_numberBlocks, M_numberBlocks);

    for (auto as : M_primalAssemblers)
    {
        unsigned int ind = as.first;
        mass.block(ind, ind).softCopy(as.second->getMass(time, sol.block(ind)));
    }

    return mass;
}

template <class InVectorType, class InMatrixType>
BlockMatrix<InMatrixType>
BlockAssembler<InVectorType, InMatrixType>::
getMassJacobian(const double& time, const BlockVector<InVectorType>& sol)
{
    BlockMatrix<InMatrixType> massJacobian;
    massJacobian.resize(M_numberBlocks, M_numberBlocks);

    for (auto as : M_primalAssemblers)
    {
        unsigned int ind = as.first;
        massJacobian.block(ind, ind).softCopy(as.second->getMassJacobian(time, sol.block(ind)));
    }

    return massJacobian;
}

template <class InVectorType, class InMatrixType>
BlockVector<InVectorType>
BlockAssembler<InVectorType, InMatrixType>::
getRightHandSide(const double& time, const BlockVector<InVectorType>& sol)
{
    BlockVector<InVectorType> rhs;
    rhs.resize(M_numberBlocks);

    for (auto as: M_primalAssemblers)
    {
        unsigned int ind = as.first;
        rhs.block(ind).softCopy(as.second->getRightHandSide(time, sol.block(ind)));
    }

    // add interface contributions
    for (auto as: M_dualAssemblers)
        as->addContributionRhs(time, rhs, sol, M_primalAssemblers.size());

    return rhs;
}

template <class InVectorType, class InMatrixType>
BlockVector<BlockVector<VectorEp>>
BlockAssembler<InVectorType, InMatrixType>::
convertFunctionRBtoFEM(BlockVector<BlockVector<DenseVector>> rbFunction,
                       EPETRACOMM comm) const
{
    BlockVector<BlockVector<VectorEp>> retVec(rbFunction.nRows());

    for (auto as : M_primalAssemblers)
    {
        unsigned int ind = as.first;
        retVec.block(ind).softCopy(as.second->convertFunctionRBtoFEM(rbFunction.block(ind)));
    }

    for (auto as : M_dualAssemblers)
    {
        unsigned int indInterface = as->getInterface().M_ID + M_primalAssemblers.size();
        retVec.block(indInterface).resize(1);
        retVec.block(indInterface).block(0).softCopy(
                VectorEp::convertDenseVector(rbFunction.block(indInterface).block(0),
                comm));
    }

    return retVec;
}

template <class InVectorType, class InMatrixType>
BlockMatrix<InMatrixType>
BlockAssembler<InVectorType, InMatrixType>::
getJacobianRightHandSide(const double& time, const BlockVector<InVectorType>& sol)
{
    BlockMatrix<InMatrixType> jac;
    jac.resize(M_numberBlocks, M_numberBlocks);

    for (auto as: M_primalAssemblers)
    {
        unsigned int ind = as.first;
        jac.block(ind,ind).softCopy(as.second->getJacobianRightHandSide(time,
                                    sol.block(ind)));
    }

    for (auto as: M_dualAssemblers)
        as->addContributionJacobianRhs(time, jac, sol, M_primalAssemblers.size());

    return jac;
}

template <class InVectorType, class InMatrixType>
std::map<unsigned int, std::string>
BlockAssembler<InVectorType, InMatrixType>::
getIDMeshTypeMap() const
{
    std::map<unsigned int, std::string> retMap;
    for (auto as: M_primalAssemblers)
    {
        std::string meshname = as.second->getTreeNode()->M_block->getMeshName();
        unsigned int dashpos = meshname.find("/");
        unsigned int formatpos = meshname.find(".mesh");
        std::string actualmeshname = meshname.substr(dashpos + 1, formatpos - dashpos - 1);
        retMap[as.first] = actualmeshname;
    }

    return retMap;
}

}
