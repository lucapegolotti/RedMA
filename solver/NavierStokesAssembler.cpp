#include <NavierStokesAssembler.hpp>

namespace RedMA
{

NavierStokesAssembler::
NavierStokesAssembler(const GetPot& datafile, const TreeNodePtr& treeNode) :
  AbstractAssembler(datafile, treeNode)
{
}

}  // namespace RedMA
