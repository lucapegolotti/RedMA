#include <AbstractAssembler.hpp>

namespace RedMA
{

AbstractAssembler::
AbstractAssembler()
{
}

AbstractAssembler::
AbstractAssembler(const TreeNodePtr& treeNode) :
  M_treeNode(treeNode)
{
}

void
AbstractAssembler::
addMapsToVector(MapVectorPtr& mapVector)
{

    for (MapVectorSTD::iterator it = M_maps.begin(); it != M_maps.end(); it++)
        mapVector->addMap(*(*it));
}

}  // namespace RedMA
