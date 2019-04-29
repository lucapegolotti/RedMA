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

    msg = std::string("Velocity FE space of size ") +
          std::to_string(M_velocityFESpace->dof().numTotalDof()) + "\n";
    printlog(GREEN, msg, M_verbose);

    std::string orderPressure = M_datafile("fluid/pressure_order", "P1");
    M_pressureFESpace.reset(new FESpace(mesh, orderPressure, 1, M_comm));

    msg = std::string("Pressure FE space of size ") +
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
    printlog(GREEN, "Assembling constant matrices ...\n", M_verbose);
    // assembleStiffnessMatrix();
    // assembleDivergenceMatrix();
    assembleMassMatrix();
    printlog(GREEN, "done\n", M_verbose);
}

void
NavierStokesAssembler::
assembleStiffnessMatrix()
{
    using namespace LifeV::ExpressionAssembly;

    printlog(YELLOW, "Assembling stiffness ...\n", M_verbose);
    M_A.reset(new Matrix(M_velocityFESpace->map()));

    double viscosity = M_datafile("fluid/viscosity", 1.0);
    integrate(elements(M_velocityFESpaceETA->mesh()),
               M_velocityFESpace->qr(),
               M_velocityFESpaceETA,
               M_velocityFESpaceETA,
               value(0.5 * viscosity) *
               dot(grad(phi_i) + transpose(grad(phi_i)),
               grad(phi_j) + transpose(grad(phi_j)))
              ) >> M_A;

    M_A->globalAssemble();
}

void
NavierStokesAssembler::
assembleDivergenceMatrix()
{
    using namespace LifeV::ExpressionAssembly;

    printlog(YELLOW, "Assembling divergence matrix ...\n", M_verbose);
    M_B.reset(new Matrix(M_pressureFESpace->map()));

    integrate(elements(M_velocityFESpaceETA->mesh()),
               M_pressureFESpace->qr(),
               M_pressureFESpaceETA,
               M_velocityFESpaceETA,
               phi_i * div(phi_j)
              ) >> M_B;

    M_B->globalAssemble(M_velocityFESpace->mapPtr(),
                        M_pressureFESpace->mapPtr());

    M_Bt.reset(new Matrix(M_velocityFESpace->map()));

    integrate(elements(M_velocityFESpaceETA->mesh()),
               M_velocityFESpace->qr(),
               M_velocityFESpaceETA,
               M_pressureFESpaceETA,
               phi_j * div(phi_i)
              ) >> M_Bt;

    M_Bt->globalAssemble(M_pressureFESpace->mapPtr(),
                         M_velocityFESpace->mapPtr());
}

void
NavierStokesAssembler::
assembleMassMatrix()
{
    using namespace LifeV::ExpressionAssembly;

    printlog(YELLOW, "Assembling mass matrix ...\n", M_verbose);
    M_M.reset(new Matrix(M_velocityFESpace->map()));

    double density = M_datafile("fluid/density", 1.0);
    integrate(elements(M_velocityFESpaceETA->mesh()),
               M_velocityFESpace->qr(),
               M_velocityFESpaceETA,
               M_velocityFESpaceETA,
               value(density) * dot (phi_i, phi_j)
              ) >> M_M;

    M_M->globalAssemble();
}

NavierStokesAssembler::MatrixPtr
NavierStokesAssembler::
getMassMatrix(const unsigned int& blockrow,
              const unsigned int& blockcol)
{
    if (blockrow == 0 && blockcol == 0)
        return M_M;

    return nullptr;
}

void
NavierStokesAssembler::
assembleConvectiveMatrix(std::vector<VectorPtr> solution)
{
    using namespace LifeV::ExpressionAssembly;

    VectorPtr velocity = solution[0];

    printlog(YELLOW, "Assembling convective matrix ...\n", M_verbose);
    M_C.reset(new Matrix(M_velocityFESpace->map()));

    double density = M_datafile("fluid/density", 1.0);
    integrate(elements(M_velocityFESpaceETA->mesh()),
               M_velocityFESpace->qr(),
               M_velocityFESpaceETA,
               M_velocityFESpaceETA,
               dot(value(density) *
                   value(M_velocityFESpaceETA , *velocity) * grad(phi_j), phi_i)
             ) >> M_C;

    M_C->globalAssemble();
}

void
NavierStokesAssembler::
assembleJacobianConvectiveMatrix(std::vector<VectorPtr> solution)
{
    using namespace LifeV::ExpressionAssembly;

    VectorPtr velocity = solution[0];

    printlog(YELLOW, "Assembling convective matrix jacobian ...\n", M_verbose);
    M_J.reset(new Matrix(M_velocityFESpace->map()));

    double density = M_datafile("fluid/density", 1.0);
    integrate(elements(M_velocityFESpaceETA->mesh()),
               M_velocityFESpace->qr(),
               M_velocityFESpaceETA,
               M_velocityFESpaceETA,
               dot(density * phi_j *
                   grad(M_velocityFESpaceETA, *velocity), phi_i)
             ) >> M_J;

    M_J->globalAssemble();
}

void
NavierStokesAssembler::
updateNonLinearTerms(const double& time, std::vector<VectorPtr> solution)
{
    assembleRhs(time);
    assembleConvectiveMatrix(solution);
    assembleJacobianConvectiveMatrix(solution);
}

void
NavierStokesAssembler::
assembleRhs(const double& time)
{
    if (!M_M)
    {
        std::string errorMsg = "Mass matrix has not been assembled yet!\n";
        throw Exception(errorMsg);
    }

    M_rhs.reset(new Vector(M_velocityFESpace->map()));
    if (!M_rhsFunction)
    {
        M_rhs->zero();
        return;
    }

    printlog(YELLOW, "Assembling right hand side ...\n", M_verbose);

    VectorPtr rhsInterp;
    rhsInterp.reset(new Vector(M_velocityFESpace->map()));

    M_velocityFESpace->interpolate(M_rhsFunction, *rhsInterp, time);
    *M_rhs = (*M_M) * (*rhsInterp);
}

NavierStokesAssembler::MatrixPtr
NavierStokesAssembler::
getJacobian(const unsigned int& blockrow, const unsigned int& blockcol)
{
    MatrixPtr retJacobian;

    if (blockrow == 0 && blockcol == 0)
    {
        retJacobian.reset(new Matrix(M_velocityFESpace->map()));
        retJacobian->zero();
        *retJacobian += *M_A;
        *retJacobian += *M_C;
        *retJacobian += *M_J;
    }
    else if (blockrow == 0 && blockcol == 1)
    {
        retJacobian = M_Bt;
    }
    else if (blockrow == 1 && blockcol == 0)
    {
        retJacobian = M_B;
    }
    else if (blockrow == 1 && blockcol == 1)
    {
        retJacobian = nullptr;
    }
    else
    {
        std::string errorMsg = "row = " + std::to_string(blockrow) + ", col = " +
                    std::to_string(blockcol) + " is an invalid combination of " +
                    "block indices for NavierStokesAssembler!";
        throw Exception(errorMsg);
    }

    return retJacobian;
}

}  // namespace RedMA
