// implementation of template class

namespace RedMA
{

template <class AssemblerType>
GlobalAssembler<AssemblerType>::
GlobalAssembler(const GetPot& datafile, commPtr_Type comm, bool verbose) :
  M_datafile(datafile),
  M_comm(comm),
  M_verbose(verbose)
{
    M_globalMap.reset(new MapEpetra());
}

template <class AssemblerType>
void
GlobalAssembler<AssemblerType>::
buildPrimalStructures(TreeStructure& tree)
{
    typedef std::map<unsigned int, TreeNodePtr>     NodesMap;

    NodesMap nodesMap = tree.getNodesMap();

    for (NodesMap::iterator it = nodesMap.begin(); it != nodesMap.end(); it++)
    {
        AssemblerTypePtr newAssembler(new AssemblerType(M_datafile, M_comm,
                                                        it->second, M_verbose));
        newAssembler->setup();

        newAssembler->addMaps(M_globalMap);
        M_assemblersMap[it->first] = newAssembler;
    }
}

template <class AssemblerType>
typename GlobalAssembler<AssemblerType>::MapEpetraPtr
GlobalAssembler<AssemblerType>::
getGlobalMap() const
{
    return M_globalMap;
}

template <class AssemblerType>
typename GlobalAssembler<AssemblerType>::MatrixPtr
GlobalAssembler<AssemblerType>::
getGlobalMass() const
{
    return M_massMatrix;
}

template <class AssemblerType>
typename GlobalAssembler<AssemblerType>::MatrixPtr
GlobalAssembler<AssemblerType>::
assembleJacobianF(const double& time, VectorPtr u) const
{
    MatrixPtr jacobian(new Matrix(*M_globalMap));

    return jacobian;
}

template <class AssemblerType>
typename GlobalAssembler<AssemblerType>::VectorPtr
GlobalAssembler<AssemblerType>::
computeF(const double& time, VectorPtr u) const
{
    VectorPtr f(new Vector(*M_globalMap));

    // implementation here
    return f;
}

template <class AssemblerType>
typename GlobalAssembler<AssemblerType>::VectorPtr
GlobalAssembler<AssemblerType>::
computeFder(const double& time, VectorPtr u) const
{
    VectorPtr fder(new Vector(*M_globalMap));

    // implementation here
    return fder;
}

template <class AssemblerType>
void
GlobalAssembler<AssemblerType>::
assembleGlobalMass()
{
    typedef std::map<unsigned int, AssemblerTypePtr> AssemblersMap;

    // M_massMatrix.reset(new Matrix(*M_globalMap));
    // M_massMatrix->zero();
    //
    // for (typename AssemblersMap::iterator it = M_assemblersMap.begin();
    //      it != M_assemblersMap.end(); it++)
    // {
    //     unsigned int blockIndex = it->first;
    //
    //     *M_massMatrix->block(blockIndex,blockIndex) =
    //         *it->second->getMassMatrix();
    // }
}

}  // namespace RedMA
