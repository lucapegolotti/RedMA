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
    M_mapVector.reset(new MapVector());
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

        newAssembler->addMapsToVector(M_mapVector);
        M_assemblersMap[it->first] = newAssembler;
    }
}

template <class AssemblerType>
typename GlobalAssembler<AssemblerType>::MapVectorPtr
GlobalAssembler<AssemblerType>::
getMapVector() const
{
    return M_mapVector;
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
    MatrixPtr jacobian(new Matrix(*M_mapVector));

    return jacobian;
}

template <class AssemblerType>
typename GlobalAssembler<AssemblerType>::VectorPtr
GlobalAssembler<AssemblerType>::
computeF(const double& time, VectorPtr u) const
{
    VectorPtr f(new Vector(*M_mapVector));

    // implementation here
    return f;
}

template <class AssemblerType>
typename GlobalAssembler<AssemblerType>::VectorPtr
GlobalAssembler<AssemblerType>::
computeFder(const double& time, VectorPtr u) const
{
    VectorPtr fder(new Vector(*M_mapVector));

    // implementation here
    return fder;
}

}  // namespace RedMA
