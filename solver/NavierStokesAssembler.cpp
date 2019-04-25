#include <NavierStokesAssembler.hpp>

namespace RedMA
{

NavierStokesAssembler::
NavierStokesAssembler() :
  AbstractAssembler()
{
}

NavierStokesAssembler::
NavierStokesAssembler(const TreeNodePtr& treeNode) :
  AbstractAssembler(treeNode)
{
}

}  // namespace RedMA
