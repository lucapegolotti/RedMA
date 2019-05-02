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
addPrimalMaps(MapEpetraPtr& globalMap, std::vector<unsigned int>& dimensions)
{
    for (MapVectorSTD::iterator it = M_primalMaps.begin();
         it != M_primalMaps.end(); it++)
    {
        *globalMap += *(*it);
        dimensions.push_back((*it)->map(LifeV::Unique)->NumGlobalElements());
    }
}

void
AbstractAssembler::
addDualMaps(MapEpetraPtr& globalMap, std::vector<unsigned int>& dimensions)
{

}

std::vector<AbstractAssembler::MapEpetraPtr>
AbstractAssembler::
getPrimalMapVector()
{
    return M_primalMaps;
}

}  // namespace RedMA
