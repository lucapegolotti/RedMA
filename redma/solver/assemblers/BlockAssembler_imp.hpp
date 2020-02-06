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
setup()
{
    typedef std::map<unsigned int, SHP(TreeNode)>         NodesMap;
    typedef aAssembler<VInner, MInner>                    InnerAssembler;
    typedef std::vector<SHP(TreeNode)>                    NodesVector;

    printlog(GREEN, "[BlockAssembler] initializing block assembler ... \n", this->M_data.getVerbose());

    NodesMap nodesMap = M_tree.getNodesMap();

    // allocate assemblers
    for (NodesMap::iterator it = nodesMap.begin(); it != nodesMap.end(); it++)
    {
        SHP(InnerAssembler) newAssembler;
        newAssembler = AssemblerFactory<VInner, MInner> (this->M_data, it->second);
        newAssembler->setup();
        M_primalAssemblers[it->second->M_ID] = newAssembler;
    }

    // allocate interface assemblers
    unsigned int interfaceID = 0;
    for (NodesMap::iterator it = nodesMap.begin(); it != nodesMap.end(); it++)
    {
        NodesVector children = it->second->M_children;

        unsigned int countOutlet = 0;
        unsigned int myID = it->second->M_ID;

        SHP(InnerAssembler) fatherAssembler = M_primalAssemblers[myID];

        unsigned int countChildren = 0;
        for (NodesVector::iterator itVector = children.begin();
             itVector != children.end(); itVector++)
        {
            if (*itVector)
            {
                unsigned int otherID = (*itVector)->M_ID;
                SHP(InnerAssembler) childAssembler = M_primalAssemblers[otherID];
                Interface<VInner, MInner> newInterface(fatherAssembler, myID,
                                                       childAssembler, otherID,
                                                       interfaceID);
                newInterface.M_indexOutlet = countChildren;
                SHP(InterfaceAssembler<VInner COMMA MInner>) inAssembler;
                inAssembler.reset(new InterfaceAssembler<VInner, MInner>(this->M_data,
                                                                         newInterface));
                M_dualAssemblers.push_back(inAssembler);
                interfaceID++;
            }
            countChildren++;
        }
    }
    M_numberBlocks = M_primalAssemblers.size() + M_dualAssemblers.size();

    printlog(GREEN, "done\n", this->M_data.getVerbose());
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
        as->addContributionRhs(rhs, sol, M_primalAssemblers.size());

    return rhs;
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
        as->addContributionJacobianRhs(jac, sol, M_primalAssemblers.size());

    return jac;
}

}
