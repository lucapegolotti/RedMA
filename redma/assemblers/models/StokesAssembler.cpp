#include "StokesAssembler.hpp"

namespace RedMA
{

StokesAssembler::
StokesAssembler(const DataContainer& data, SHP(TreeNode) treeNode) :
  M_dataContainer(data),
  M_treeNode(treeNode)
{
    M_density = data("fluid/density", 1.0);
    M_viscosity = data("fluid/viscosity", 0.035);
    M_name = "StokesAssembler";
}


SHP(BlockVector)
StokesAssembler::
buildZeroVector() const
{
    SHP(VECTOREPETRA) uComp(new VECTOREPETRA(M_velocityFESpace->map(),
                                             LifeV::Unique));
    uComp->zero();

    SHP(VECTOREPETRA) pComp(new VECTOREPETRA(M_pressureFESpace->map(),
                                             LifeV::Unique));

    pComp->zero();

    SHP(BlockVector) retVec;
    retVec->resize(2);
    retVec->block(0)->setData(uComp);
    retVec->block(1)->setData(pComp);
    return retVec;
}

SHP(BlockVector)
StokesAssembler::
getForcingTerm(const double& time) const
{
    // we don't consider forcing terms for the moment
    return buildZeroVector();
}

void
StokesAssembler::
addNeumannBCs(SHP(BlockVector)& input, const double& time,
              const SHP(BlockVector)& sol)
{
    // // if (this->M_treeNode->isOutletNode())
    //     // //     this->M_bcManager->applyNeumannBc(time, input, M_velocityFESpace, 0, flowRates);
    //     //
    //     // if (this->M_treeNode->isOutletNode())
    //     // {
    //     //     auto flowRates = computeFlowRates(sol, true);
    //     //
    //     //     for (auto rate : flowRates)
    //     //     {
    //     //         double P = this->M_bcManager->getNeumannBc(time, rate.first, rate.second);
    //     //         BlockVector<InVectorType> curvec(this->M_nComponents);
    //     //
    //     //         curvec.block(0).data().reset(new VECTOREPETRA(*M_flowRateVectors[rate.first]));
    //     //         curvec.block(0) *= P;
    //     //         input += curvec;
    //     //
    //     //         // addBackFlowStabilization(input, sol, rate.first);
    //     //     }
    //     // }
}

SHP(BlockMatrix)
StokesAssembler::
assembleReducedStiffness(SHP(BCManager) bcManager,
                         BlockMDEIMStructure* structure)
{
    using namespace LifeV;
    using namespace ExpressionAssembly;

    SHP(BlockMatrix) stiffness(new BlockMatrix());

    stiffness->resize(2, 2);
    bool useFullStrain = M_dataContainer("fluid/use_strain", true);

    SHP(MATRIXEPETRA) A(new MATRIXEPETRA(M_velocityFESpace->map()));

    if (structure)
    {
        unsigned int numVolumes = (*structure)(0,0)->numReducedElements;
        unsigned int* volumes = (*structure)(0,0)->reducedElements.data();

        if (useFullStrain)
        {
            integrate(elements(M_velocityFESpaceETA->mesh(), 0, numVolumes, volumes, true),
                      M_velocityFESpace->qr(),
                      M_velocityFESpaceETA,
                      M_velocityFESpaceETA,
                      value(0.5 * M_viscosity) *
                      dot(grad(phi_i) + transpose(grad(phi_i)),
                      grad(phi_j) + transpose(grad(phi_j)))
                  ) >> A;
        }
        else
        {
            integrate(elements(M_velocityFESpaceETA->mesh(), 0, numVolumes, volumes, true),
                      M_velocityFESpace->qr(),
                      M_velocityFESpaceETA,
                      M_velocityFESpaceETA,
                      value(M_viscosity) *
                      dot(grad(phi_i),grad(phi_j))
                  ) >> A;
        }
    }
    else
    {
        if (useFullStrain)
        {
            integrate(elements(M_velocityFESpaceETA->mesh()),
                      M_velocityFESpace->qr(),
                      M_velocityFESpaceETA,
                      M_velocityFESpaceETA,
                      value(0.5 * M_viscosity) *
                      dot(grad(phi_i) + transpose(grad(phi_i)),
                      grad(phi_j) + transpose(grad(phi_j)))
                  ) >> A;
        }
        else
        {
            integrate(elements(M_velocityFESpaceETA->mesh()),
                      M_velocityFESpace->qr(),
                      M_velocityFESpaceETA,
                      M_velocityFESpaceETA,
                      value(M_viscosity) *
                      dot(grad(phi_i),grad(phi_j))
                  ) >> A;
        }
    }
    A->globalAssemble();

    stiffness->block(0,0)->setData(A);

    bcManager->apply0DirichletMatrix(*stiffness, M_velocityFESpace, 0, 0.0);

    return stiffness;
}

SHP(BlockMatrix)
StokesAssembler::
assembleReducedMass(SHP(BCManager) bcManager,
                    BlockMDEIMStructure* structure)
{
    using namespace LifeV;
    using namespace ExpressionAssembly;

    SHP(BlockMatrix) mass(new BlockMatrix());

    mass->resize(2, 2);

    SHP(MATRIXEPETRA) M(new MATRIXEPETRA(M_velocityFESpace->map()));

    if (structure)
    {
        unsigned int numVolumes = (*structure)(0,0)->numReducedElements;
        unsigned int* volumes = (*structure)(0,0)->reducedElements.data();
        integrate(elements(M_velocityFESpaceETA->mesh(), 0, numVolumes, volumes, true),
                  M_velocityFESpace->qr(),
                  M_velocityFESpaceETA,
                  M_velocityFESpaceETA,
                  value(M_density) * dot(phi_i, phi_j)
              ) >> M;

    }
    else
    {
        integrate(elements(M_velocityFESpaceETA->mesh()),
                  M_velocityFESpace->qr(),
                  M_velocityFESpaceETA,
                  M_velocityFESpaceETA,
                  value(M_density) * dot(phi_i, phi_j)
              ) >> M;
    }
    M->globalAssemble();
    mass->block(0,0)->setData(M);

    bcManager->apply0DirichletMatrix(*mass, M_velocityFESpace, 0, 1.0);

    return mass;
}

SHP(BlockMatrix)
StokesAssembler::
assembleReducedDivergence(SHP(BCManager) bcManager,
                          BlockMDEIMStructure* structure)
{
    using namespace LifeV;
    using namespace ExpressionAssembly;

    SHP(BlockMatrix) divergence(new BlockMatrix());

    divergence->resize(2, 2);

    SHP(MATRIXEPETRA) BT(new MATRIXEPETRA(this->M_velocityFESpace->map()));

    if (structure)
    {
        unsigned int numVolumes = (*structure)(0,1)->numReducedElements;
        unsigned int* volumes = (*structure)(0,1)->reducedElements.data();
        integrate(elements(M_velocityFESpaceETA->mesh(), 0, numVolumes, volumes, true),
                  M_velocityFESpace->qr(),
                  M_velocityFESpaceETA,
                  M_pressureFESpaceETA,
                  value(-1.0) * phi_j * div(phi_i)
              ) >> BT;
    }
    else
    {
        integrate(elements(M_velocityFESpaceETA->mesh()),
                  M_velocityFESpace->qr(),
                  M_velocityFESpaceETA,
                  M_pressureFESpaceETA,
                  value(-1.0) * phi_j * div(phi_i)
              ) >> BT;
    }

    BT->globalAssemble(M_pressureFESpace->mapPtr(),
                       M_velocityFESpace->mapPtr());

    SHP(MATRIXEPETRA) B(new MATRIXEPETRA(M_pressureFESpace->map()));

    if (structure)
    {
        unsigned int numVolumes = (*structure)(1,0)->numReducedElements;
        unsigned int* volumes = (*structure)(1,0)->reducedElements.data();
        integrate(elements(M_velocityFESpaceETA->mesh(), 0, numVolumes, volumes, true),
                 M_pressureFESpace->qr(),
                 M_pressureFESpaceETA,
                 M_velocityFESpaceETA,
                 phi_i * div(phi_j)
             ) >> B;
    }
    else
    {
        integrate(elements(M_velocityFESpaceETA->mesh()),
                 M_pressureFESpace->qr(),
                 M_pressureFESpaceETA,
                 M_velocityFESpaceETA,
                 phi_i * div(phi_j)
             ) >> B;
    }
    B->globalAssemble(M_velocityFESpace->mapPtr(),
                      M_pressureFESpace->mapPtr());

    divergence->block(0,1)->setData(BT);
    divergence->block(1,0)->setData(B);

    bcManager->apply0DirichletMatrix(*divergence, M_velocityFESpace, 0, 0.0);

    return divergence;
}

void
StokesAssembler::
assembleFlowRateVectors()
{
    // assemble inflow flow rate vector
    if (M_treeNode->isInletNode())
    {
        auto face = M_treeNode->M_block->getInlet();
        // std::cout << "\nInlet normal\n" << std::endl;
        // face.print();
        M_flowRateVectors[face.M_flag] = assembleFlowRateVector(face);
    }

    if (M_treeNode->isOutletNode())
    {
        auto faces = M_treeNode->M_block->getOutlets();

        for (auto face : faces)
            M_flowRateVectors[face.M_flag] = assembleFlowRateVector(face);
    }
}

void
StokesAssembler::
assembleFlowRateJacobians(SHP(BCManager) bcManager)
{
    // assemble inflow flow rate vector
    if (M_treeNode->isInletNode())
    {
        auto face = M_treeNode->M_block->getInlet();
        M_flowRateJacobians[face.M_flag]->resize(2,2);
        M_flowRateJacobians[face.M_flag]->block(0,0)->setData(assembleFlowRateJacobian(face));

        apply0DirichletBCsMatrix(bcManager,M_flowRateJacobians[face.M_flag], 0.0);
    }

    if (M_treeNode->isOutletNode())
    {
        auto faces = M_treeNode->M_block->getOutlets();

        for (auto face : faces)
        {
            M_flowRateJacobians[face.M_flag]->resize(2,2);
            M_flowRateJacobians[face.M_flag]->block(0,0)->setData(assembleFlowRateJacobian(face));

            apply0DirichletBCsMatrix(bcManager,M_flowRateJacobians[face.M_flag], 0.0);
        }
    }
}

void
StokesAssembler::
apply0DirichletBCsMatrix(SHP(BCManager) bcManager,
                         SHP(BlockMatrix)& matrix, double diagCoeff)
{
    bcManager->apply0DirichletMatrix(*matrix, M_velocityFESpace,
                                     0, diagCoeff);
}

std::map<unsigned int, double>
StokesAssembler::
computeFlowRates(SHP(BlockVector) sol, bool verbose)
{
    std::string msg;
    std::map<unsigned int, double> flowRates;
    if (M_treeNode->isInletNode())
    {
        auto face = M_treeNode->M_block->getInlet();

        flowRates[face.M_flag] = static_cast<VECTOREPETRA*>(sol->block(0)->data().get())->dot(*M_flowRateVectors[face.M_flag]);
        std::string msg = "[";
        msg += this->M_name;
        msg += "]  inflow rate = ";
        msg += std::to_string(flowRates[face.M_flag]);
        msg += "\n";
        printlog(YELLOW, msg, verbose);
    }

    if (this->M_treeNode->isOutletNode())
    {
        auto faces = this->M_treeNode->M_block->getOutlets();

        for (auto face : faces)
        {
            flowRates[face.M_flag] = static_cast<VECTOREPETRA*>(sol->block(0)->data().get())->dot(*M_flowRateVectors[face.M_flag]);
            std::string msg = "[";
            msg += this->M_name;
            msg += "]  outflow rate = ";
            msg += std::to_string(flowRates[face.M_flag]);
            msg += "\n";
            printlog(YELLOW, msg, verbose);
        }
    }

    return flowRates;
}

SHP(VECTOREPETRA)
StokesAssembler::
assembleFlowRateVector(const GeometricFace& face)
{
    using namespace LifeV;
    using namespace ExpressionAssembly;

    SHP(VECTOREPETRA) flowRateVectorRepeated;
    flowRateVectorRepeated.reset(new VECTOREPETRA(M_velocityFESpace->map(),
                                                  Repeated));

    QuadratureBoundary myBDQR(buildTetraBDQR(quadRuleTria7pt));

    integrate(boundary(M_velocityFESpaceETA->mesh(), face.M_flag),
              myBDQR,
              M_velocityFESpaceETA,
              dot(phi_i, Nface)
          ) >> flowRateVectorRepeated;

    flowRateVectorRepeated->globalAssemble();

    SHP(VECTOREPETRA) flowRateVector(new VECTOREPETRA(*flowRateVectorRepeated,
                                                      Unique));
    return flowRateVector;
}

void
StokesAssembler::
addBackFlowStabilization(SHP(BlockVector)& input,
                         SHP(BlockVector) sol,
                         const unsigned int& faceFlag)
{
    // using namespace LifeV;
    // using namespace ExpressionAssembly;
    //
    // SHP(VECTOREPETRA) vn(new VECTOREPETRA(*sol.block(0).data()));
    //
    // *vn *= *M_flowRateVectors[faceFlag];
    //
    // SHP(VECTOREPETRA) absvn(new VECTOREPETRA(*vn));
    // absvn->abs();
    //
    // *vn -= *absvn;
    // *vn /= 2.0;
    //
    // *vn *= *sol.block(0).data();
    //
    // SHP(VECTOREPETRA) vnRepeated(new VECTOREPETRA(*vn, Repeated));
    // SHP(VECTOREPETRA) backflowStabRepeated(new VECTOREPETRA(vn->mapPtr(), Repeated));
    //
    // QuadratureBoundary myBDQR(buildTetraBDQR(quadRuleTria7pt));
    //
    // integrate(boundary(M_velocityFESpaceETA->mesh(), faceFlag),
    //           myBDQR,
    //           M_velocityFESpaceETA,
    //           dot(value(M_velocityFESpaceETA, *vnRepeated), phi_i)
    //       ) >> backflowStabRepeated;
    //
    // backflowStabRepeated->globalAssemble();
    //
    // *backflowStabRepeated *= 0.2 * M_density;
    //
    // SHP(VECTOREPETRA) backflowStab(new VECTOREPETRA(*backflowStabRepeated, Unique));
    //
    // *input.block(0).data() += *backflowStab;
}

void
StokesAssembler::
exportNorms(double t, SHP(VECTOREPETRA) velocity, SHP(VECTOREPETRA) pressure)
{
    bool exportNorms = M_dataContainer("exporter/exportnorms", true);

    if (exportNorms)
    {
        std::string outputName = M_dataContainer("exporter/outdir", "solutions/") + "block";
        outputName += std::to_string(this->M_treeNode->M_ID);
        outputName += "_norm.txt";
        std::ofstream filename(outputName, std::ios_base::app);
        filename << t << ",";
        filename << M_velocityFESpace->h1Norm(*velocity) << ",";
        filename << M_pressureFESpace->l2Norm(*pressure) << "\n";
    }
}

void
StokesAssembler::
initializePythonStructures()
{
    // setenv("PYTHONPATH",".",1);

    // Py_Initialize();
    // PyObject* pName = PyUnicode_DecodeFSDefault("test");

    // M_pModule = PyImport_Import(pName);
    // Py_DECREF(pName);

    // M_pFunc = PyObject_GetAttrString(M_pModule, "evaluate_model");
}

void
StokesAssembler::
computeWallShearStress(SHP(VECTOREPETRA) velocity, SHP(VECTOREPETRA) WSS,
                       EPETRACOMM comm)
{
    using namespace LifeV;
    using namespace ExpressionAssembly;

    QuadratureBoundary myBDQR(buildTetraBDQR(quadRuleTria7pt));

    unsigned int wallFlag = 10;
    if (M_massWall == nullptr)
    {
        M_massWall.reset(new MATRIXEPETRA(M_velocityFESpace->map()));

        integrate(boundary(M_velocityFESpaceETA->mesh(), wallFlag),
                  myBDQR,
                  M_velocityFESpaceETA,
                  M_velocityFESpaceETA,
                  dot(phi_i, phi_j)
              ) >> M_massWall;
        M_massWall->globalAssemble();
    }

    WSS->zero();

    SHP(VECTOREPETRA) velocityRepeated(new VECTOREPETRA(*velocity,
                                                         Repeated));
    SHP(VECTOREPETRA) weakWSSRepeated(new VECTOREPETRA(M_velocityFESpace->map(), Repeated));

    integrate(boundary(M_velocityFESpaceETA->mesh(), wallFlag),
              myBDQR,
              M_velocityFESpaceETA,
              value(M_viscosity) *
              dot(
             //  (dot(grad(M_velocityFESpaceETA, *velocityRepeated) +
             //  transpose(grad(M_velocityFESpaceETA, *velocityRepeated)),Nface)
             // //  -
             // // dot(dot(grad(M_velocityFESpaceETA, *velocityRepeated) +
             // //  transpose(grad(M_velocityFESpaceETA, *velocityRepeated))
             // //  ,Nface),Nface) * Nface
             // )
             (grad(M_velocityFESpaceETA, *velocityRepeated) +
             transpose(grad(M_velocityFESpaceETA, *velocityRepeated))) * Nface
             -
             dot(
             (grad(M_velocityFESpaceETA, *velocityRepeated) +
             transpose(grad(M_velocityFESpaceETA, *velocityRepeated))) * Nface,
             Nface
             ) * Nface,
             phi_i)
             ) >> weakWSSRepeated;

    weakWSSRepeated->globalAssemble();
    SHP(VECTOREPETRA) weakWSSUnique(new VECTOREPETRA(*weakWSSRepeated, Unique));

    LinearSolver linearSolver(comm);
    linearSolver.setOperator(M_massWall);

    Teuchos::RCP<Teuchos::ParameterList> aztecList =
                                    Teuchos::rcp(new Teuchos::ParameterList);
    aztecList = Teuchos::getParametersFromXmlFile("datafiles/SolverParamList.xml");
    linearSolver.setParameters(*aztecList);

    typedef LifeV::PreconditionerML         precML_type;
    typedef std::shared_ptr<precML_type>    precMLPtr_type;
    precML_type * precRawPtr;
    precRawPtr = new precML_type;

    GetPot dummyDatafile;
    precRawPtr->setDataFromGetPot(dummyDatafile, "precMLL");
    std::shared_ptr<LifeV::Preconditioner> precPtr;
    precPtr.reset(precRawPtr);

    linearSolver.setPreconditioner(precPtr);
    linearSolver.setRightHandSide(weakWSSUnique);
    linearSolver.solve(WSS);
}

void
StokesAssembler::
initializeVelocityFESpace(EPETRACOMM comm)
{
    // initialize fespace velocity
    std::string orderVelocity = M_dataContainer("fluid/velocity_order", "P2");
    M_velocityFESpace.reset(new FESPACE(this->M_treeNode->M_block->getMesh(),
                                       orderVelocity, 3, comm));

    M_velocityFESpaceETA.reset(new ETFESPACE3(M_velocityFESpace->mesh(),
                                           &(M_velocityFESpace->refFE()),
                                             comm));
}


void
StokesAssembler::
initializePressureFESpace(EPETRACOMM comm)
{
    // initialize fespace pressure
    std::string orderPressure = M_dataContainer("fluid/pressure_order", "P1");

    M_pressureFESpace.reset(new FESPACE(this->M_treeNode->M_block->getMesh(),
                                       orderPressure, 1, comm));
    M_pressureFESpaceETA.reset(new ETFESPACE1(M_pressureFESpace->mesh(),
                                           &(M_pressureFESpace->refFE()),
                                             comm));
}

}
