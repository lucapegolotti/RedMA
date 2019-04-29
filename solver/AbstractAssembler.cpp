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
addMaps(MapEpetraPtr& globalMap, std::vector<unsigned int>& dimensions)
{
    for (MapVectorSTD::iterator it = M_maps.begin(); it != M_maps.end(); it++)
    {
        *globalMap += *(*it);
        dimensions.push_back((*it)->map(LifeV::Unique)->NumGlobalElements());
    }
}

std::vector<AbstractAssembler::MapEpetraPtr>
AbstractAssembler::
getMapVector()
{
    return M_maps;
}

}  // namespace RedMA
