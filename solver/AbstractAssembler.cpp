#include <AbstractAssembler.hpp>

namespace RedMA
{

AbstractAssembler::
AbstractAssembler(const GetPot& datafile, commPtr_Type comm,
                  const TreeNodePtr& treeNode, bool verbose) :
  M_datafile(datafile),
  M_comm(comm),
  M_treeNode(treeNode),
  M_verbose(verbose)
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
