#include <AbstractAssembler.hpp>

namespace RedMA
{

AbstractAssembler::
AbstractAssembler(const GetPot& datafile, const TreeNodePtr& treeNode) :
  M_datafile(datafile),
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
