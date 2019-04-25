#include <NavierStokesAssembler.hpp>

namespace RedMA
{

NavierStokesAssembler::
NavierStokesAssembler(const GetPot& datafile, commPtr_Type comm,
                      const TreeNodePtr& treeNode) :
  AbstractAssembler(datafile, comm, treeNode)
{
}

void
NavierStokesAssembler::
setup()
{
    MeshPtr mesh = M_treeNode->M_block->getMesh();
    std::string orderVelocity = M_datafile("fluid/velocity_order", "P2");
    M_velocityFESpace.reset(new FESpace(mesh, orderVelocity, 3, M_comm));
}

}  // namespace RedMA
