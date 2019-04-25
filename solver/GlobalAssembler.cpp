#include <GlobalAssembler.hpp>

namespace RedMA
{

template <class AssemblerType>
GlobalAssembler<AssemblerType>::
GlobalAssembler()
{

}

template <class AssemblerType>
void
GlobalAssembler<AssemblerType>::
buildPrimalStructures(TreeStructure& tree, MapVectorPtr& mapVector,
                      MatrixStructuredPtr& matrixStructured)
{
    typedef std::map<unsigned int, TreeNodePtr>     NodesMap;

    NodesMap nodesMap = tree.getNodesMap();

    for (NodesMap::iterator it = nodesMap.begin(); it != nodesMap.end(); it++)
    {
        AssemblerTypePtr newAssembler(new AssemblerType(it->second));
        newAssembler->addMapsToVector(mapVector);
        M_assemblersMap[it->first] = newAssembler;
    }
}

}  // namespace RedMA
