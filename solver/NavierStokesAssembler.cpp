#include <NavierStokesAssembler.hpp>

namespace RedMA
{

NavierStokesAssembler::
NavierStokesAssembler(const GetPot& datafile, commPtr_Type comm,
                      const TreeNodePtr& treeNode, bool verbose) :
  AbstractAssembler(datafile, comm, treeNode, verbose)
{
}

void
NavierStokesAssembler::
setup()
{
    std::string msg = std::string("[NavierStokesAssembler] setting up assembler");
    msg = msg + " with building block of type " + M_treeNode->M_block->name() +
          " ...\n";
    printlog(MAGENTA, msg, M_verbose);
    MeshPtr mesh = M_treeNode->M_block->getMesh();

    std::string orderVelocity = M_datafile("fluid/velocity_order", "P2");
    M_velocityFESpace.reset(new FESpace(mesh, orderVelocity, 3, M_comm));

    msg = std::string("[NavierStokesAssembler] velocity FE space of size ") +
          std::to_string(M_velocityFESpace->dof().numTotalDof()) + "\n";
    printlog(GREEN, msg, M_verbose);

    std::string orderPressure = M_datafile("fluid/velocity_pressure", "P1");
    M_pressureFESpace.reset(new FESpace(mesh, orderPressure, 1, M_comm));

    msg = std::string("[NavierStokesAssembler] pressure FE space of size ") +
          std::to_string(M_pressureFESpace->dof().numTotalDof()) + "\n";
    printlog(GREEN, msg, M_verbose);

    // add maps to vector
    M_maps.push_back(M_velocityFESpace->mapPtr());
    M_maps.push_back(M_pressureFESpace->mapPtr());

    // fespaces
    M_velocityFESpaceETA.reset(new ETFESpaceVelocity(M_velocityFESpace->mesh(),
                                                   &(M_velocityFESpace->refFE()),
                                                     M_comm));

    M_pressureFESpaceETA.reset(new ETFESpacePressure(M_pressureFESpace->mesh(),
                                                   &(M_pressureFESpace->refFE()),
                                                     M_comm));

    assembleConstantMatrices();
    printlog(MAGENTA, "done\n", M_verbose);
}

void
NavierStokesAssembler::
assembleConstantMatrices()
{
    assembleStiffness();
}

void
NavierStokesAssembler::
assembleStiffness()
{
    M_stiffness.reset(new Matrix(M_velocityFESpace->map()));

    double viscosity = M_datafile("fluid/viscosity", 1.0);
    using namespace LifeV::ExpressionAssembly;

    integrate(elements(M_velocityFESpaceETA->mesh()),
              M_velocityFESpace->qr(),
              M_velocityFESpaceETA,
              M_velocityFESpaceETA,
              value(0.5 * viscosity) *
              dot(grad(phi_i) + transpose(grad(phi_i)),
              grad(phi_j) + transpose(grad(phi_j)))
    ) >> M_stiffness;

    M_stiffness->globalAssemble();
}

void
NavierStokesAssembler::
assembleDivergenceMatrix()
{
    using namespace LifeV::ExpressionAssembly;

}

void
NavierStokesAssembler::
assembleMass()
{
    using namespace LifeV::ExpressionAssembly;

}

}  // namespace RedMA
