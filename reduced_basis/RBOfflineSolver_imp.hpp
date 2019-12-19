namespace RedMA
{

template<class AssemblerType>
RBOfflineSolver<AssemblerType>::
RBOfflineSolver(const GetPot& datafile, rbLifeV::ParameterHandler& parameterHandler,
                commPtr_Type comm, bool verbose) :
  rbLifeV::ApproximatedAffineFemProblem(comm, parameterHandler),
  GlobalSolver<AssemblerType>(datafile, comm, verbose)
{
}

template<class AssemblerType>
const std::shared_ptr<LifeV::MapEpetra>
RBOfflineSolver<AssemblerType>::
getFieldMap(int iField) const
{
    return this->M_globalAssembler.getLocalMap(iField);
}

template<class AssemblerType>
typename RBOfflineSolver<AssemblerType>::FESpacePtr
RBOfflineSolver<AssemblerType>::
getFieldFeSpace(int iField) const
{
    return this->M_globalAssembler.getLocalFespace(iField);
}

template<class AssemblerType>
void
RBOfflineSolver<AssemblerType>::
setParameter(const param_Type& mu)
{
    if (mu != M_myCurrentMu)
    {
        this->M_tree.resetMeshes();
        unsigned int offset = 0;
        setParameterSubdomains(mu, offset);

        this->M_tree.traverseAndDeformGeometries();
        this->M_globalAssembler.setPhysicalParameters(mu, offset);

        // TODO: maybe there is something more to be done here

    }
}

template<class AssemblerType>
void
RBOfflineSolver<AssemblerType>::
setGeometricalParameterSubdomains(param_Type mu, unsigned int& offset)
{
    typedef std::shared_ptr<TreeNode>   TreeNodePtr;
    typedef std::vector<TreeNodePtr>    TreeNodesVector;

    TreeNodePtr root = this->M_tree.getRoot();

    std::queue<TreeNodePtr> nodesQueue;
    nodesQueue.push(root);
    offset = 0;
    while(nodesQueue.size() != 0)
    {
        TreeNodePtr curNode = nodesQueue.front();
        nodesQueue.pop();

        curNode->M_block->setGeometricalParameters(mu, offset);

        TreeNodesVector& children = curNode->M_children;
        unsigned int expectedChildren =
                     curNode->M_block->expectedNumberOfChildren();
        for (int i = 0; i < expectedChildren; i++)
        {
            TreeNodePtr curChild = children[i];

            if (curChild)
            {
                GeometricFace curFace = curNode->M_block->getOutlet(i);
                nodesQueue.push(curChild);
            }
        }
    }

}

const
typedef RBOfflineSolver<AssemblerType>::MeshPtr
RBOfflineSolver<AssemblerType>::
mesh() const
{
    new throw Exception("mesh() method not implemented!");
    return nullptr;
}

}
