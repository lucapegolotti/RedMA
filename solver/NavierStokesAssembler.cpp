#include <NavierStokesAssembler.hpp>

namespace RedMA
{

NavierStokesAssembler::
NavierStokesAssembler(const GetPot& datafile, commPtr_Type comm,
                      const TreeNodePtr& treeNode, bool verbose) :
  AbstractAssembler(datafile, comm, treeNode, verbose)
{
    M_addNoslipBC = true;
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

    M_couplingFESpace.reset(new FESpace(mesh, orderVelocity, 3, M_comm));

    msg = std::string("Velocity FE space of size ") +
          std::to_string(M_velocityFESpace->dof().numTotalDof() * 3) + "\n";
    printlog(GREEN, msg, M_verbose);

    std::string orderPressure = M_datafile("fluid/pressure_order", "P1");
    M_pressureFESpace.reset(new FESpace(mesh, orderPressure, 1, M_comm));

    msg = std::string("Pressure FE space of size ") +
          std::to_string(M_pressureFESpace->dof().numTotalDof()) + "\n";
    printlog(GREEN, msg, M_verbose);

    // add maps to vector
    M_primalMaps.push_back(M_velocityFESpace->mapPtr());
    M_primalMaps.push_back(M_pressureFESpace->mapPtr());

    // fespaces
    M_velocityFESpaceETA.reset(new ETFESpaceVelocity(M_velocityFESpace->mesh(),
                                                   &(M_velocityFESpace->refFE()),
                                                     M_comm));

    M_pressureFESpaceETA.reset(new ETFESpacePressure(M_pressureFESpace->mesh(),
                                                   &(M_pressureFESpace->refFE()),
                                                     M_comm));

    M_couplingFESpaceETA.reset(new ETFESpaceCoupling(M_velocityFESpace->mesh(),
                                                   &(M_velocityFESpace->refFE()),
                                                     M_comm));

    M_indexCoupling = 0;

    assembleConstantMatrices();
    setExporter();

    M_useStabilization = M_datafile("fluid/use_stabilization", false);
    printlog(MAGENTA, "done\n", M_verbose);
}

void
NavierStokesAssembler::
assembleConstantMatrices()
{
    printlog(GREEN, "Assembling constant matrices ...\n", M_verbose);
    assembleStiffnessMatrix();
    assembleDivergenceMatrix();
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
    bool useFullStrain = M_datafile("fluid/use_strain", false);

    if (useFullStrain)
        integrate(elements(M_velocityFESpaceETA->mesh()),
                   M_velocityFESpace->qr(),
                   M_velocityFESpaceETA,
                   M_velocityFESpaceETA,
                   value(0.5 * viscosity) *
                   dot(grad(phi_i) + transpose(grad(phi_i)),
                   grad(phi_j) + transpose(grad(phi_j)))
                  ) >> M_A;
    else
        integrate(elements(M_velocityFESpaceETA->mesh()),
                 M_velocityFESpace->qr(),
                 M_velocityFESpaceETA,
                 M_velocityFESpaceETA,
                 value(viscosity) *
                 dot(grad(phi_i),grad(phi_j))
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
               value(-1.0) * phi_j * div(phi_i)
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
    {
        // we copy the matrix so that we cannot change M_M (e.g. by applying
        // boundary conditions)
        MatrixPtr returnMatrix(new Matrix(*M_M));
        return returnMatrix;
    }

    return nullptr;
}

void
NavierStokesAssembler::
assembleConvectiveMatrix()
{
    using namespace LifeV::ExpressionAssembly;

    VectorPtr velocityRepeated(new Vector(*M_prevSolution[0], LifeV::Repeated));

    printlog(YELLOW, "Assembling convective matrix ...\n", M_verbose);
    if (M_C == nullptr)
        M_C.reset(new Matrix(M_velocityFESpace->map()));
    else
        M_C->zero();

    double density = M_datafile("fluid/density", 1.0);
    integrate(elements(M_velocityFESpaceETA->mesh()),
               M_velocityFESpace->qr(),
               M_velocityFESpaceETA,
               M_velocityFESpaceETA,
               dot(value(density) *
                   value(M_velocityFESpaceETA , *velocityRepeated) * grad(phi_j),
                         phi_i)
             ) >> M_C;

    M_C->globalAssemble();
}

void
NavierStokesAssembler::
assembleJacobianConvectiveMatrix()
{
    using namespace LifeV::ExpressionAssembly;

    VectorPtr velocityRepeated(new Vector(*M_prevSolution[0], LifeV::Repeated));

    printlog(YELLOW, "Assembling convective matrix jacobian ...\n", M_verbose);
    if (M_J == nullptr)
        M_J.reset(new Matrix(M_velocityFESpace->map()));
    else
        M_J->zero();

    double density = M_datafile("fluid/density", 1.0);
    integrate(elements(M_velocityFESpaceETA->mesh()),
               M_velocityFESpace->qr(),
               M_velocityFESpaceETA,
               M_velocityFESpaceETA,
               dot(density * phi_j *
                   grad(M_velocityFESpaceETA, *velocityRepeated), phi_i)
             ) >> M_J;

    M_J->globalAssemble();
}

void
NavierStokesAssembler::
setTimeAndPrevSolution(const double& time, std::vector<VectorPtr> solution,
                       bool assembleBlocks)
{
    M_time = time;
    M_prevSolution = solution;
    // reset all things that need to be recomputed

    if (assembleBlocks)
    {
        assembleConvectiveMatrix();
        assembleJacobianConvectiveMatrix();
        assembleForcingTerm();
        assembleForcingTermTimeDerivative();
    }

    if (M_useStabilization)
    {
        if (M_stabilization == nullptr)
        {
            std::string velocityOrder = M_datafile("fluid/velocity_order", "P1");
            unsigned int uOrder;
            if (!std::strcmp(velocityOrder.c_str(),"P1"))
                uOrder = 1;
            else if (!std::strcmp(velocityOrder.c_str(),"P2"))
                uOrder = 2;
            else
                throw Exception("Implement suitable velocity order");

            M_stabilization.reset(new
                         VMS_SUPGStabilization(this->M_timeIntegrationOrder,
                                               uOrder,
                                               M_velocityFESpace,
                                               M_pressureFESpace,
                                               M_velocityFESpaceETA,
                                               M_pressureFESpaceETA));
            double density = M_datafile("fluid/density", 1.0);
            double viscosity = M_datafile("fluid/viscosity", 1.0);
            M_stabilization->setDensityAndViscosity(density, viscosity);
        }
        if (assembleBlocks)
            M_stabilization->assembleBlocks(M_prevSolution[0], M_prevSolution[1],
                                            M_forcingTerm, this->M_dt);
    }
}

void
NavierStokesAssembler::
assembleForcingTerm()
{
    if (!M_M)
    {
        std::string errorMsg = "Mass matrix has not been assembled yet!\n";
        throw Exception(errorMsg);
    }

    M_forcingTerm.reset(new Vector(M_velocityFESpace->map()));
    if (!M_forceFunction)
    {
        M_forcingTerm->zero();
        return;
    }

    printlog(YELLOW, "Assembling right hand side ...\n", M_verbose);

    VectorPtr forcingTermInterp;
    forcingTermInterp.reset(new Vector(M_velocityFESpace->map()));

    M_velocityFESpace->interpolate(M_forceFunction, *forcingTermInterp, M_time);
    *M_forcingTerm = (*M_M) * (*forcingTermInterp);
}

NavierStokesAssembler::MatrixPtr
NavierStokesAssembler::
getJacobian(const unsigned int& blockrow, const unsigned int& blockcol)
{
    std::string blockStr = std::string("(") + std::to_string(blockrow) + "," +
                           std::to_string(blockcol) + ")";

    printlog(GREEN, "Constructing jacobian, block " + blockStr +
                    " ...\n", M_verbose);

    MatrixPtr retJacobian;

    if (blockrow == 0 && blockcol == 0)
    {
        retJacobian.reset(new Matrix(M_velocityFESpace->map()));
        retJacobian->zero();
        *retJacobian += *M_A;
        *retJacobian += *M_C;
        *retJacobian += *M_J;
        if (M_useStabilization)
            *retJacobian += (*M_stabilization->block00Jac());
        retJacobian->globalAssemble();
    }
    else if (blockrow == 0 && blockcol == 1)
    {
        retJacobian.reset(new Matrix(M_velocityFESpace->map()));
        *retJacobian = *M_Bt;

        if (M_useStabilization)
            *retJacobian += (*M_stabilization->block01Jac());

        retJacobian->globalAssemble(M_pressureFESpace->mapPtr(),
                                    M_velocityFESpace->mapPtr());
    }
    else if (blockrow == 1 && blockcol == 0)
    {
        retJacobian.reset(new Matrix(M_pressureFESpace->map()));
        *retJacobian = *M_B;

        if (M_useStabilization)
            *retJacobian += (*M_stabilization->block10Jac());

        retJacobian->globalAssemble(M_velocityFESpace->mapPtr(),
                                    M_pressureFESpace->mapPtr());
    }
    else if (blockrow == 1 && blockcol == 1)
    {
        retJacobian = nullptr;
        if (M_useStabilization)
        {
            retJacobian.reset(new Matrix(*M_stabilization->block11Jac()));
            retJacobian->globalAssemble();
        }
    }
    else if (blockrow > numberOfBlocks() || blockcol > numberOfBlocks())
    {
        std::string errorMsg = "row = " + std::to_string(blockrow) + ", col = " +
                    std::to_string(blockcol) + " is an invalid combination of " +
                    "block indices for NavierStokesAssembler!";
        throw Exception(errorMsg);
    }
    // we multiply by -1 because it is at the right hand side when writing
    // H du/dt = F(t,u)
    if (retJacobian)
        *retJacobian *= (-1.0);
    return retJacobian;
}

std::vector<NavierStokesAssembler::VectorPtr>
NavierStokesAssembler::
computeF()
{
    std::vector<VectorPtr> Fs;
    VectorPtr velocity = M_prevSolution[0];
    VectorPtr pressure = M_prevSolution[1];

    // assemble F first component
    VectorPtr F1;
    F1.reset(new Vector(M_velocityFESpace->map()));
    F1->zero();

    *F1 += *M_forcingTerm;
    *F1 -= (*M_A) * (*velocity);
    *F1 -= (*M_C) * (*velocity);
    *F1 -= (*M_Bt) * (*pressure);

    if (M_useStabilization)
    {
        *F1 -= *M_stabilization->velocityResidual(velocity, pressure,
                                                  M_forcingTerm, this->M_dt);
    }

    unsigned int count = numberOfBlocks();
    for (std::map<unsigned int, MatrixPtr>::iterator it = getMapsQTs().begin();
         it != getMapsQTs().end(); it++)
    {
        *F1 -= (*it->second) * (*M_prevSolution[count]);
        count++;
    }

    applyNeumannBCs(F1, M_inflowLaw);

    // assemble F second component
    VectorPtr F2;
    F2.reset(new Vector(M_pressureFESpace->map()));
    F2->zero();

    *F2 = (*M_B) * (*velocity);
    *F2 *= (-1);

    if (M_useStabilization)
    {
        *F2 -= *M_stabilization->pressureResidual(velocity, pressure,
                                                  M_forcingTerm, this->M_dt);
    }

    Fs.push_back(F1);
    Fs.push_back(F2);

    count = 0;
    for (std::map<unsigned int, MatrixPtr>::iterator it = getMapsQs().begin();
         it != getMapsQs().end(); it++)
    {
        VectorPtr newF;
        MatrixPtr curCouplingMatrix = it->second;
        newF.reset(new Vector(*M_dualMaps[count]));
        newF->zero();

        *newF -= (*curCouplingMatrix) * (*velocity);
        Fs.push_back(newF);
        count++;
    }
    return Fs;
}

void
NavierStokesAssembler::
assembleForcingTermTimeDerivative()
{
    if (!M_M)
    {
        std::string errorMsg = "Mass matrix has not been assembled yet!\n";
        throw Exception(errorMsg);
    }

    M_forcingTermTimeDer.reset(new Vector(M_velocityFESpace->map()));
    if (!M_forceTimeDerFunction)
    {
        M_forcingTermTimeDer->zero();
        return;
    }

    printlog(YELLOW, "Assembling right hand side derivative ...\n", M_verbose);

    VectorPtr forcingTermDerInterp;
    forcingTermDerInterp.reset(new Vector(M_velocityFESpace->map()));

    M_velocityFESpace->interpolate(M_forceTimeDerFunction, *forcingTermDerInterp,
                                   M_time);
    *M_forcingTermTimeDer = (*M_M) * (*forcingTermDerInterp);
}

std::vector<NavierStokesAssembler::VectorPtr>
NavierStokesAssembler::
computeFder()
{
    std::vector<VectorPtr> Fs;
    VectorPtr velocity = M_prevSolution[0];
    VectorPtr pressure = M_prevSolution[1];

    // assemble F first component
    VectorPtr F1;
    F1.reset(new Vector(M_velocityFESpace->map()));
    F1->zero();

    *F1 += *M_forcingTermTimeDer;

    applyNeumannBCs(F1, M_inflowLawDt);

    // assemble F second component
    VectorPtr F2;
    F2.reset(new Vector(M_pressureFESpace->map()));
    F2->zero();

    Fs.push_back(F1);
    Fs.push_back(F2);
    unsigned int count = 0;
    for (std::map<unsigned int, MatrixPtr>::iterator it = getMapsQs().begin();
         it != getMapsQs().end(); it++)
    {
        VectorPtr newF;
        MatrixPtr curCouplingMatrix = it->second;
        newF.reset(new Vector(*M_dualMaps[count]));
        newF->zero();
        Fs.push_back(newF);
        count++;
    }
    return Fs;
}

NavierStokesAssembler::MatrixPtr
NavierStokesAssembler::
getUpdateMass(const unsigned int& blockrow, const unsigned int& blockcol)
{
    if (M_useStabilization)
    {
        if (blockrow == 0 && blockcol == 0)
            return M_stabilization->blockMass00();
        if (blockrow == 1 && blockcol == 0)
            return M_stabilization->blockMass10();
        if (blockrow == 0 && blockcol == 1)
            return M_stabilization->blockMass01();
        if (blockrow == 1 && blockcol == 1)
            return M_stabilization->blockMass11();
    }
    return nullptr;
}

NavierStokesAssembler::MatrixPtr
NavierStokesAssembler::
getUpdateMassJac(const unsigned int& blockrow, const unsigned int& blockcol)
{
    if (M_useStabilization)
    {
        if (blockrow == 0 && blockcol == 0)
            return M_stabilization->blockMass00Jac();
        if (blockrow == 1 && blockcol == 0)
            return M_stabilization->blockMass10Jac();
        if (blockrow == 0 && blockcol == 1)
            return M_stabilization->blockMass01Jac();
        if (blockrow == 1 && blockcol == 1)
            return M_stabilization->blockMass11Jac();
    }
    return nullptr;
}

NavierStokesAssembler::MatrixPtr
NavierStokesAssembler::
getUpdateMassJacVelocity(const unsigned int& blockrow, const unsigned int& blockcol)
{
    if (M_useStabilization)
    {
        if (blockrow == 0 && blockcol == 0)
            return M_stabilization->assembleMassWithVelocity(M_prevSolution[0],
                                                             this->M_dt);
    }
    return nullptr;
}

NavierStokesAssembler::BoundaryConditionPtr
NavierStokesAssembler::
createBCHandler(std::function<double(double)> law)
{
    LifeV::BCFunctionBase zeroFunction (fZero);

    BoundaryConditionPtr bcs;
    bcs.reset(new LifeV::BCHandler);

    const unsigned int inletFlag = 1;
    const unsigned int wallFlag = 10;
    const unsigned int inletRing = 30;
    const unsigned int outletRing = 31;

    // if this is the root node, we impose dirichlet boundary conditions at
    // inlet too (if inflow conditions are dirichlet)
    std::string inflowBCType = M_datafile("fluid/inflow_bc","dirichlet");
    if (M_treeNode->M_ID == 0)
    {
        if (std::strcmp(inflowBCType.c_str(),"dirichlet") == 0)
        {

            auto inflowBoundaryCondition = std::bind(poiseulleInflow,
                                                     std::placeholders::_1,
                                                     std::placeholders::_2,
                                                     std::placeholders::_3,
                                                     std::placeholders::_4,
                                                     std::placeholders::_5,
                                                     M_treeNode->M_block->getInlet(),
                                                     law);

            LifeV::BCFunctionBase inflowFunction(inflowBoundaryCondition);
            bcs->addBC("Inlet", inletFlag, LifeV::Essential, LifeV::Full,
                        inflowFunction, 3);
        }
    }

    if (M_addNoslipBC)
        bcs->addBC("Wall", wallFlag, LifeV::Essential,
                    LifeV::Full, zeroFunction, 3);
    // clamping the root and leaf nodes
    else
    {
        if (M_treeNode->M_ID == 0)
            bcs->addBC("InletRing", inletRing, LifeV::EssentialEdges,
                       LifeV::Full, zeroFunction,   3);
        if (M_treeNode->M_nChildren == 0)
            bcs->addBC("OutletRing", outletRing, LifeV::EssentialEdges,
                       LifeV::Full, zeroFunction,   3);
    }
    return bcs;
}

double
NavierStokesAssembler::
fZero(const double& t, const double& x, const double& y, const double& z,
      const unsigned int& i)
{
    return 0.0;
}

double
NavierStokesAssembler::
fOne(const double& t, const double& x, const double& y, const double& z,
      const unsigned int& i)
{
    return 1.0;
}

double
NavierStokesAssembler::
poiseulleInflow(const double& t, const double& x, const double& y,
                const double& z, const unsigned int& i,
                const GeometricFace& face,
                std::function<double(double)> maxLaw)
{
    typedef LifeV::VectorSmall<3>   Vector3D;

    const Vector3D& center = face.M_center;
    const Vector3D& normal = face.M_normal;
    double radius = face.M_radius;

    Vector3D curPoint(x,y,z);
    Vector3D diff = curPoint - center;
    double normDiff = radius - diff.norm();

    double inflowNorm = maxLaw(t) * normDiff * normDiff / (radius * radius);
    Vector3D inflow = -inflowNorm * normal;
    return inflow[i];
}

double
NavierStokesAssembler::
neumannInflow(const double& t, const double& x, const double& y,
              const double& z, const unsigned int& i,
              std::function<double(double)> maxLaw)
{
    return maxLaw(t);
}

void
NavierStokesAssembler::
setLawInflow(std::function<double(double)> maxLaw)
{
    M_inflowLaw = maxLaw;
}

void
NavierStokesAssembler::
setLawDtInflow(std::function<double(double)> maxLawDt)
{
    M_inflowLawDt = maxLawDt;
}

void
NavierStokesAssembler::
applyNeumannBCs(VectorPtr vector, std::function<double(double)> law)
{
    BoundaryConditionPtr bcs;
    bcs.reset(new LifeV::BCHandler);

    const unsigned int inletFlag = 1;
    const unsigned int wallFlag = 10;

    std::string inflowBCType = M_datafile("fluid/inflow_bc","dirichlet");
    if (M_treeNode->M_ID == 0 &&
        std::strcmp(inflowBCType.c_str(),"neumann") == 0)
    {

        auto inflowBoundaryCondition = std::bind(neumannInflow,
                                                 std::placeholders::_1,
                                                 std::placeholders::_2,
                                                 std::placeholders::_3,
                                                 std::placeholders::_4,
                                                 std::placeholders::_5,
                                                 law);

        LifeV::BCFunctionBase inflowFunction(inflowBoundaryCondition);
        bcs->addBC("Inlet", inletFlag, LifeV::Natural, LifeV::Normal,
                    inflowFunction);
    }

    updateBCs(bcs, M_velocityFESpace);

    bcManageRhs(*vector, *M_velocityFESpace->mesh(), M_velocityFESpace->dof(),
                *bcs, M_velocityFESpace->feBd(), 1.0, M_time);
}

void
NavierStokesAssembler::
applyBCsRhsRosenbrock(std::vector<VectorPtr> rhs,
                      std::vector<VectorPtr> utilde,
                      const double& time,
                      const double& dt,
                      const double& alphai,
                      const double& gammai)
{
    BoundaryConditionPtr bc = createBCHandler(M_inflowLaw);

    updateBCs(bc, M_velocityFESpace);

    VectorPtr auxVec(new Vector(M_velocityFESpace->map()));
    auxVec->zero();

    bcManageRhs(*auxVec, *M_velocityFESpace->mesh(), M_velocityFESpace->dof(),
                *bc, M_velocityFESpace->feBd(), 1.0, time + dt * alphai);

    BoundaryConditionPtr bcDt = createBCHandler(M_inflowLawDt);

    updateBCs(bcDt, M_velocityFESpace);

    VectorPtr auxVecDt(new Vector(M_velocityFESpace->map()));
    auxVecDt->zero();

    bcManageRhs(*auxVecDt, *M_velocityFESpace->mesh(), M_velocityFESpace->dof(),
                *bcDt, M_velocityFESpace->feBd(), 1.0, time);

    VectorPtr rhsValues(new Vector(M_velocityFESpace->map()));
    rhsValues->zero();

    *rhsValues += *auxVecDt;
    *rhsValues *= (gammai * dt);
    *rhsValues += *auxVec;
    *rhsValues -= *utilde[0];

    BoundaryConditionPtr finalBcs(new LifeV::BCHandler);
    LifeV::BCVector bcVectorDirichlet(*rhsValues,
                                      M_velocityFESpace->map().mapSize()/3);

    const unsigned int inletFlag = 1;
    std::string inflowBCType = M_datafile("fluid/inflow_bc","dirichlet");
    if (M_treeNode->M_ID == 0 &&
        std::strcmp(inflowBCType.c_str(),"dirichlet") == 0)
    {
        finalBcs->addBC("Inflow", inletFlag, LifeV::Essential,
                        LifeV::Full, bcVectorDirichlet, 3);
    }

    LifeV::BCFunctionBase zeroFunction (fZero);

    const unsigned int wallFlag = 10;

    finalBcs->addBC("Wall", wallFlag, LifeV::Essential,
                    LifeV::Full, zeroFunction, 3);

    updateBCs(finalBcs, M_velocityFESpace);
    bcManageRhs(*rhs[0], *M_velocityFESpace->mesh(), M_velocityFESpace->dof(),
                *finalBcs, M_velocityFESpace->feBd(), 1.0, time);
}

void
NavierStokesAssembler::
applyBCsBackwardEuler(std::vector<VectorPtr> rhs, const double& coeff,
                      const double& time)
{
    BoundaryConditionPtr bc = createBCHandler(M_inflowLaw);
    updateBCs(bc, M_velocityFESpace);

    bcManageRhs(*rhs[0], *M_velocityFESpace->mesh(), M_velocityFESpace->dof(),
                *bc, M_velocityFESpace->feBd(), coeff, time);
}

void
NavierStokesAssembler::
updateBCs(BoundaryConditionPtr bcToUpdate, FESpacePtr fespace)
{
    bcToUpdate->bcUpdate(*fespace->mesh(), fespace->feBd(), fespace->dof());
}

void
NavierStokesAssembler::
applyBCsMatrix(MatrixPtr matrix, const double& diagonalCoefficient,
               const unsigned int& iblock, const unsigned int& jblock)
{
    std::string msg = "[NavierStokesAssembler] applying boundary conditions to";
    msg += " block (" + std::to_string(iblock) + "," + std::to_string(jblock) +
            ") ...\n";
    printlog(MAGENTA, msg, M_verbose);

    if (matrix)
    {
        MapEpetra rangeMap = matrix->rangeMap();
        MapEpetra domainMap = matrix->domainMap();

        BoundaryConditionPtr bc = createBCHandler(M_inflowLaw);
        updateBCs(bc, M_velocityFESpace);

        if (iblock == 0 && jblock == 0)
        {
            bcManageMatrix(*matrix, *M_velocityFESpace->mesh(),
                           M_velocityFESpace->dof(), *bc, M_velocityFESpace->feBd(),
                           diagonalCoefficient, 0.0);
        }
        if (iblock == 0 && jblock == 1)
        {
            bcManageMatrix(*matrix, *M_velocityFESpace->mesh(),
                           M_velocityFESpace->dof(), *bc, M_velocityFESpace->feBd(),
                           0.0, 0.0);
        }
        matrix->globalAssemble(std::make_shared<MapEpetra>(domainMap),
                               std::make_shared<MapEpetra>(rangeMap));
    }
}

void
NavierStokesAssembler::
setExporter()
{
    std::string outputName = "block";
    outputName += std::to_string(M_treeNode->M_ID);

    std::string outdir = M_datafile("exporter/outdir", "solutions/");
    boost::filesystem::create_directory(outdir);

    std::string format = M_datafile("exporter/type", "hdf5");
    if (!std::strcmp(format.c_str(), "hdf5"))
        M_exporter.reset(new ExporterHDF5(M_datafile, outputName));
    else
        M_exporter.reset(new ExporterVTK(M_datafile, outputName));
    M_exporter->setMeshProcId(M_velocityFESpace->mesh(), M_comm->MyPID());

    M_velocityExporter.reset(new Vector(M_velocityFESpace->map(),
                                        M_exporter->mapType()));
    M_pressureExporter.reset(new Vector(M_pressureFESpace->map(),
                                        M_exporter->mapType()));
    M_lagrangeMultiplierExporter.reset(new Vector(M_velocityFESpace->map(),
                                       M_exporter->mapType()));

    M_exporter->addVariable(LifeV::ExporterData<Mesh>::VectorField,
                         "velocity", M_velocityFESpace, M_velocityExporter, 0.0);

    M_exporter->addVariable(LifeV::ExporterData<Mesh>::ScalarField,
                         "pressure", M_pressureFESpace, M_pressureExporter, 0.0);

    M_exporter->addVariable(LifeV::ExporterData<Mesh>::VectorField,
                         "lagrange_multiplier", M_velocityFESpace, M_lagrangeMultiplierExporter, 0.0);

    M_exporter->setPostDir(outdir);
}

void
NavierStokesAssembler::
exportSolutions(const double& time, std::vector<VectorPtr> solutions)
{
    printlog(MAGENTA, "[NavierStokesAssembler] exporting solution ...\n",
             M_verbose);
    *M_velocityExporter = *solutions[0];
    *M_pressureExporter = *solutions[1];
    *M_lagrangeMultiplierExporter = *reconstructLagrangeMultipliers(solutions, 2);
    CoutRedirecter ct;
    ct.redirect();
    M_exporter->postProcess(time);
    printlog(CYAN, ct.restore(), M_verbose);
    printlog(MAGENTA, "done\n", M_verbose);
}

std::vector<double>
NavierStokesAssembler::
computeNorms(std::vector<VectorPtr> solutions)
{
    printlog(MAGENTA, "[NavierStokesAssembler] writing norm solution ...\n",
             M_verbose);

    std::vector<double> norms;

    VectorPtr velocityRepeated(new Vector(*solutions[0], LifeV::Repeated));
    norms.push_back(M_velocityFESpace->l2Norm(*velocityRepeated));
    norms.push_back(M_velocityFESpace->h1Norm(*velocityRepeated));
    VectorPtr pressureRepeated(new Vector(*solutions[1], LifeV::Repeated));
    norms.push_back(M_pressureFESpace->l2Norm(*pressureRepeated));

    printlog(MAGENTA, "done\n", M_verbose);
    return norms;
}

void
NavierStokesAssembler::
checkResidual(std::vector<VectorPtr> solutions,
              std::vector<VectorPtr> prevSolutions, double dt)
{
    std::cout << "Checking residuals" << std::endl;

    VectorPtr F1;
    F1.reset(new Vector(M_velocityFESpace->map()));
    F1->zero();

    *F1 = *M_forcingTerm;
    *F1 *= (-1);
    *F1 = (*M_A) * (*solutions[0]);

    M_prevSolution[0] = solutions[0];

    assembleConvectiveMatrix();

    *F1 += (*M_C) * (*solutions[0]);
    *F1 += (*M_Bt) * (*solutions[1]);
    *F1 *= (dt);

    VectorPtr Fmass;
    Fmass.reset(new Vector(M_velocityFESpace->map()));
    Fmass->zero();

    *Fmass = (*M_M) * (*solutions[0]);
    *Fmass -= (*M_M) * (*prevSolutions[0]);

    *F1 += *Fmass;

    std::vector<VectorPtr> ff;
    ff.push_back(F1);
    applyBCsBackwardEuler(ff,0,0);
    std::cout << "Residual 1 = " << F1->norm2() << std::endl << std::flush;

    if (M_useStabilization)
    {
        M_stabilization->assembleBlocks(solutions[0], solutions[1],
                                        M_forcingTerm, dt);
        *F1 += (*M_stabilization->blockMass00()) * (*solutions[0]);
        *F1 -= (*M_stabilization->blockMass00()) * (*prevSolutions[0]);
        *F1 += (dt) * (*M_stabilization->velocityResidual(solutions[0], solutions[1],
                                             M_forcingTerm, dt));

        applyBCsBackwardEuler(ff,0,0);
        std::cout << "Residual 1 after stabilization = " << F1->norm2() << std::endl << std::flush;
    }
    else
    {
        std::string velocityOrder = M_datafile("fluid/velocity_order", "P1");
        unsigned int uOrder;
        if (!std::strcmp(velocityOrder.c_str(),"P1"))
            uOrder = 1;
        else if (!std::strcmp(velocityOrder.c_str(),"P2"))
            uOrder = 2;
        else
            throw Exception("Implement suitable velocity order");

        M_stabilization.reset(new
                     VMS_SUPGStabilization(this->M_timeIntegrationOrder,
                                           uOrder,
                                           M_velocityFESpace,
                                           M_pressureFESpace,
                                           M_velocityFESpaceETA,
                                           M_pressureFESpaceETA));
        double density = M_datafile("fluid/density", 1.0);
        double viscosity = M_datafile("fluid/viscosity", 1.0);
        M_stabilization->setDensityAndViscosity(density, viscosity);

        M_stabilization->assembleBlocks(solutions[0], solutions[1],
                                        M_forcingTerm, this->M_dt);

        VectorPtr fstab;
        fstab.reset(new Vector(M_velocityFESpace->map()));
        fstab->zero();

        *fstab += (*M_stabilization->blockMass00()) * (*solutions[0]);
        *fstab -= (*M_stabilization->blockMass00()) * (*prevSolutions[0]);
        *fstab += (dt) * (*M_stabilization->velocityResidual(solutions[0], solutions[1],
                                             M_forcingTerm, dt));

        std::cout << "Residual of the stabilization 1 = " << fstab->norm2() << std::endl;
    }
    VectorPtr F2;
    F2.reset(new Vector(M_pressureFESpace->map()));
    F2->zero();

    *F2 += (*M_B) * (*solutions[0]);
    *F2 *= (dt);

    std::cout << "Residual 2 = " << F2->norm2() << std::endl << std::flush;

    if (M_useStabilization)
    {
        *F2 += (*M_stabilization->blockMass10()) * (*solutions[0]);
        *F2 -= (*M_stabilization->blockMass10()) * (*prevSolutions[0]);
        *F2 += (dt) * (*M_stabilization->pressureResidual(solutions[0], solutions[1],
                                             M_forcingTerm, dt));

        std::cout << "Residual 2 after stabilization = " << F2->norm2() << std::endl << std::flush;
    }
    else
    {
        VectorPtr fstab;
        fstab.reset(new Vector(M_pressureFESpace->map()));
        fstab->zero();

        *fstab += (*M_stabilization->blockMass10()) * (*solutions[0]);
        *fstab -= (*M_stabilization->blockMass10()) * (*prevSolutions[0]);
        *fstab += (dt) * (*M_stabilization->pressureResidual(solutions[0], solutions[1],
                                             M_forcingTerm, dt));

        std::cout << "Residual of the stabilization 2 = " << fstab->norm2() << std::endl;
    }
    std::cout << "Total residual = " << std::sqrt(F1->norm2() * F1->norm2() +  F2->norm2() * F2->norm2()) << std::endl << std::flush;
}

std::string
NavierStokesAssembler::
normFileFirstLine()
{
    return "time,ul2,uh1,l2p\n";
}


}  // namespace RedMA
