#include <PseudoFSI.hpp>

namespace RedMA
{

PseudoFSI::
PseudoFSI(const GetPot& datafile, commPtr_Type comm,
          const TreeNodePtr& treeNode, bool verbose) :
  NavierStokesAssembler(datafile, comm, treeNode, verbose)
{
}

void
PseudoFSI::
assembleMassMatrix()
{
    NavierStokesAssembler::assembleMassMatrix();
}

}  // namespace RedMA