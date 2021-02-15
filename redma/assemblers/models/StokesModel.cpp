#include "StokesModel.hpp"

namespace RedMA
{

StokesModel::
StokesModel(const DataContainer& data, shp<TreeNode> treeNode) :
  M_dataContainer(data),
  M_treeNode(treeNode)
{
    M_density = data("fluid/density", 1.0);
    M_viscosity = data("fluid/viscosity", 0.035);
    M_velocityOrder = data("fluid/velocity_order", "P2");
    M_pressureOrder = data("fluid/pressure_order", "P1");
    M_addNoSlipBC = true;
}


shp<BlockVector>
StokesModel::
buildZeroVector() const
{
    shp<VECTOREPETRA> uComp(new VECTOREPETRA(M_velocityFESpace->map(),
                                             LifeV::Unique));
    uComp->zero();
    shp<VECTOREPETRA> pComp(new VECTOREPETRA(M_pressureFESpace->map(),
                                             LifeV::Unique));

    pComp->zero();

    shp<BlockVector> retVec(new BlockVector(2));

    retVec->setBlock(0,wrap(uComp));
    retVec->setBlock(1,wrap(pComp));

    return retVec;
}

shp<aVector>
StokesModel::
getForcingTerm(const double& time) const
{
    // we don't consider forcing terms for the moment
    return buildZeroVector();
}

shp<aMatrix>
StokesModel::
assembleReducedStiffness(shp<BCManager> bcManager)
{
    using namespace LifeV;
    using namespace ExpressionAssembly;

    shp<BlockMatrix> stiffness(new BlockMatrix(2,2));

    bool useFullStrain = M_dataContainer("fluid/use_strain", true);

    shp<MATRIXEPETRA> A(new MATRIXEPETRA(M_velocityFESpace->map()));

    // if (structure)
    if (0)
    {
        // unsigned int numVolumes = (*structure)(0,0)->numReducedElements;
        // unsigned int* volumes = (*structure)(0,0)->reducedElements.data();
        //
        // if (useFullStrain)
        // {
        //     integrate(elements(M_velocityFESpaceETA->mesh(), 0, numVolumes, volumes, true),
        //               M_velocityFESpace->qr(),
        //               M_velocityFESpaceETA,
        //               M_velocityFESpaceETA,
        //               value(0.5 * M_viscosity) *
        //               dot(grad(phi_i) + transpose(grad(phi_i)),
        //               grad(phi_j) + transpose(grad(phi_j)))
        //           ) >> A;
        // }
        // else
        // {
        //     integrate(elements(M_velocityFESpaceETA->mesh(), 0, numVolumes, volumes, true),
        //               M_velocityFESpace->qr(),
        //               M_velocityFESpaceETA,
        //               M_velocityFESpaceETA,
        //               value(M_viscosity) *
        //               dot(grad(phi_i),grad(phi_j))
        //           ) >> A;
        // }
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

    shp<SparseMatrix> Awrapper(new SparseMatrix);
    Awrapper->setData(A);
    stiffness->setBlock(0,0,Awrapper);

    bcManager->apply0DirichletMatrix(*stiffness, M_velocityFESpace, 0, 0.0, !(M_addNoSlipBC));

    return stiffness;
}

shp<aMatrix>
StokesModel::
assembleReducedMass(shp<BCManager> bcManager)
{
    using namespace LifeV;
    using namespace ExpressionAssembly;

    shp<BlockMatrix> mass(new BlockMatrix(2,2));
    shp<MATRIXEPETRA> M(new MATRIXEPETRA(M_velocityFESpace->map()));

    // if (structure)
    if (0)
    {
        // unsigned int numVolumes = (*structure)(0,0)->numReducedElements;
        // unsigned int* volumes = (*structure)(0,0)->reducedElements.data();
        // integrate(elements(M_velocityFESpaceETA->mesh(), 0, numVolumes, volumes, true),
        //           M_velocityFESpace->qr(),
        //           M_velocityFESpaceETA,
        //           M_velocityFESpaceETA,
        //           value(M_density) * dot(phi_i, phi_j)
        //       ) >> M;

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
    shp<SparseMatrix> Mwrapper(new SparseMatrix);
    Mwrapper->setData(M);
    mass->setBlock(0,0,Mwrapper);

    bcManager->apply0DirichletMatrix(*mass, M_velocityFESpace, 0, 1.0, !(M_addNoSlipBC));

    return mass;
}

shp<aMatrix>
StokesModel::
assembleReducedDivergence(shp<BCManager> bcManager)
{
    using namespace LifeV;
    using namespace ExpressionAssembly;

    shp<BlockMatrix> divergence(new BlockMatrix(2,2));

    shp<MATRIXEPETRA> BT(new MATRIXEPETRA(this->M_velocityFESpace->map()));

    // if (structure)
    if (false)
    {
        // unsigned int numVolumes = (*structure)(0,1)->numReducedElements;
        // unsigned int* volumes = (*structure)(0,1)->reducedElements.data();
        // integrate(elements(M_velocityFESpaceETA->mesh(), 0, numVolumes, volumes, true),
        //           M_velocityFESpace->qr(),
        //           M_velocityFESpaceETA,
        //           M_pressureFESpaceETA,
        //           value(-1.0) * phi_j * div(phi_i)
        //       ) >> BT;
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

    shp<MATRIXEPETRA> B(new MATRIXEPETRA(M_pressureFESpace->map()));

    // if (structure)
    if (false)
    {
        // unsigned int numVolumes = (*structure)(1,0)->numReducedElements;
        // unsigned int* volumes = (*structure)(1,0)->reducedElements.data();
        // integrate(elements(M_velocityFESpaceETA->mesh(), 0, numVolumes, volumes, true),
        //          M_pressureFESpace->qr(),
        //          M_pressureFESpaceETA,
        //          M_velocityFESpaceETA,
        //          phi_i * div(phi_j)
        //      ) >> B;
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

    shp<SparseMatrix> BTwrapper(new SparseMatrix);
    BTwrapper->setData(BT);

    shp<SparseMatrix> Bwrapper(new SparseMatrix);
    Bwrapper->setData(B);

    divergence->setBlock(0,1, BTwrapper);
    divergence->setBlock(1,0, Bwrapper);

    bcManager->apply0DirichletMatrix(*divergence, M_velocityFESpace, 0, 0.0, !(M_addNoSlipBC));

    return divergence;
}

void
StokesModel::
assembleFlowRateVectors()
{
    // assemble inflow flow rate vector
    if (M_treeNode->isInletNode())
    {
        auto face = M_treeNode->M_block->getInlet();

        M_flowRateVectors[face.M_flag] = assembleFlowRateVector(face);
    }

    // assemble outflow(s) flow rate vector(s)
    if (M_treeNode->isOutletNode())
    {
        auto faces = M_treeNode->M_block->getOutlets();

        for (auto face : faces)
            M_flowRateVectors[face.M_flag] = assembleFlowRateVector(face);
    }
}

void
StokesModel::
assembleFlowRateJacobians(shp<BCManager> bcManager)
{
    // assemble inflow flow rate vector
    // if (M_treeNode->isInletNode())
    // {
    //     auto face = M_treeNode->M_block->getInlet();
    //     M_flowRateJacobians[face.M_flag]->resize(2,2);
    //     M_flowRateJacobians[face.M_flag]->block(0,0)->setData(assembleFlowRateJacobian(face));
    //
    //     applyDirichletBCsMatrix(bcManager,M_flowRateJacobians[face.M_flag], 0.0);
    // }

    if (M_treeNode->isOutletNode())
    {
        auto faces = M_treeNode->M_block->getOutlets();

        for (auto face : faces)
        {
            shp<BlockMatrix> newJacobian(new BlockMatrix(2,2));
            shp<SparseMatrix> jacWrapper(new SparseMatrix());
            jacWrapper->setMatrix(assembleFlowRateJacobian(face));

            newJacobian->setBlock(0,0,jacWrapper);

            applyDirichletBCsMatrix(bcManager, newJacobian, 0.0);

            M_flowRateJacobians[face.M_flag] = newJacobian;

        }
    }
}

void
StokesModel::
applyDirichletBCsMatrix(shp<BCManager> bcManager,
                         shp<aMatrix> matrix, double diagCoeff)
{
    auto matrixConverted = spcast<BlockMatrix>(matrix);
    bcManager->apply0DirichletMatrix(*matrixConverted, M_velocityFESpace,
                                     0, diagCoeff, !(M_addNoSlipBC));
}

std::map<unsigned int, double>
StokesModel::
computeFlowRates(shp<aVector> sol, bool verbose)
{
    auto solBlck = convert<BlockVector>(sol);

    std::string msg;
    std::map<unsigned int, double> flowRates;
    if (M_treeNode->isInletNode())
    {
        auto face = M_treeNode->M_block->getInlet();

        flowRates[face.M_flag] = static_cast<VECTOREPETRA*>(solBlck->block(0)->data().get())->dot(*M_flowRateVectors[face.M_flag]);
        std::string msg = "[";
        msg += "StokesModel";
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
            flowRates[face.M_flag] = static_cast<VECTOREPETRA*>(solBlck->block(0)->data().get())->dot(*M_flowRateVectors[face.M_flag]);
            std::string msg = "[";
            msg += "StokesModel";
            msg += "]  outflow rate = ";
            msg += std::to_string(flowRates[face.M_flag]);
            msg += "\n";
            printlog(YELLOW, msg, verbose);
        }
    }

    return flowRates;
}

shp<VECTOREPETRA>
StokesModel::
assembleFlowRateVector(const GeometricFace& face)
{
    using namespace LifeV;
    using namespace ExpressionAssembly;

    shp<VECTOREPETRA> flowRateVectorRepeated;
    flowRateVectorRepeated.reset(new VECTOREPETRA(M_velocityFESpace->map(),
                                                  Repeated));

    QuadratureBoundary myBDQR(buildTetraBDQR(quadRuleTria7pt));

    integrate(boundary(M_velocityFESpaceETA->mesh(), face.M_flag),
              myBDQR,
              M_velocityFESpaceETA,
              dot(phi_i, Nface)
          ) >> flowRateVectorRepeated;

    flowRateVectorRepeated->globalAssemble();

    shp<VECTOREPETRA> flowRateVector(new VECTOREPETRA(*flowRateVectorRepeated,
                                                      Unique));
    return flowRateVector;
}

shp<MATRIXEPETRA>
StokesModel::
assembleFlowRateJacobian(const GeometricFace& face)
{
    using namespace LifeV;
    using namespace ExpressionAssembly;

    const double dropTolerance(2.0 * std::numeric_limits<double>::min());

    shp<MAPEPETRA> rangeMap = M_flowRateVectors[face.M_flag]->mapPtr();
    EPETRACOMM comm = rangeMap->commPtr();

    Epetra_Map epetraMap = M_flowRateVectors[face.M_flag]->epetraMap();
    unsigned int numElements = epetraMap.NumMyElements();
    unsigned int numGlobalElements = epetraMap.NumGlobalElements();

    // TODO: this should be optimized
    shp<MATRIXEPETRA> flowRateJacobian;
    flowRateJacobian.reset(new MATRIXEPETRA(M_velocityFESpace->map(), numGlobalElements, false));

    // compute outer product of flowrate vector with itself
    for (unsigned int j = 0; j < numGlobalElements; j++)
    {
        double myvaluecol = 0;

        if (M_flowRateVectors[face.M_flag]->isGlobalIDPresent(j))
            myvaluecol = M_flowRateVectors[face.M_flag]->operator[](j);

        double valuecol = 0;
        comm->SumAll(&myvaluecol, &valuecol, 1);

        if (std::abs(valuecol) > dropTolerance)
        {
            for (unsigned int i = 0; i < numElements; i++)
            {
                unsigned int gdof = epetraMap.GID(i);
                if (M_flowRateVectors[face.M_flag]->isGlobalIDPresent(gdof))
                {
                    double valuerow = M_flowRateVectors[face.M_flag]->operator[](gdof);
                    if (std::abs(valuerow * valuecol) > dropTolerance)
                    {
                        flowRateJacobian->addToCoefficient(gdof, j, valuerow * valuecol);
                    }
                }
            }
        }

    }

    comm->Barrier();

    flowRateJacobian->globalAssemble();

    return flowRateJacobian;
}

void
StokesModel::
addBackFlowStabilization(shp<aVector>& input,
                         shp<aVector> sol,
                         const unsigned int& faceFlag)
{
    // using namespace LifeV;
    // using namespace ExpressionAssembly;
    //
    // shp<VECTOREPETRA> vn(new VECTOREPETRA(*sol.block(0).data()));
    //
    // *vn *= *M_flowRateVectors[faceFlag];
    //
    // shp<VECTOREPETRA> absvn(new VECTOREPETRA(*vn));
    // absvn->abs();
    //
    // *vn -= *absvn;
    // *vn /= 2.0;
    //
    // *vn *= *sol.block(0).data();
    //
    // shp<VECTOREPETRA> vnRepeated(new VECTOREPETRA(*vn, Repeated));
    // shp<VECTOREPETRA> backflowStabRepeated(new VECTOREPETRA(vn->mapPtr(), Repeated));
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
    // shp<VECTOREPETRA> backflowStab(new VECTOREPETRA(*backflowStabRepeated, Unique));
    //
    // *input.block(0).data() += *backflowStab;
}

void
StokesModel::
exportNorms(const double& t, shp<VECTOREPETRA> velocity, shp<VECTOREPETRA> pressure)
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
StokesModel::
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
StokesModel::
computeWallShearStress(shp<VECTOREPETRA> velocity, shp<VECTOREPETRA> WSS,
                       EPETRACOMM comm)
{
    using namespace LifeV;
    using namespace ExpressionAssembly;

    QuadratureBoundary myBDQR(buildTetraBDQR(quadRuleTria7pt));

    unsigned int wallFlag = M_treeNode->M_block->wallFlag();
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

    shp<VECTOREPETRA> velocityRepeated(new VECTOREPETRA(*velocity,
                                                         Repeated));
    shp<VECTOREPETRA> weakWSSRepeated(new VECTOREPETRA(M_velocityFESpace->map(), Repeated));

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
    shp<VECTOREPETRA> weakWSSUnique(new VECTOREPETRA(*weakWSSRepeated, Unique));

    LinearSolver linearSolver(comm);
    linearSolver.setOperator(M_massWall);

    Teuchos::RCP<Teuchos::ParameterList> aztecList =
                                    Teuchos::rcp(new Teuchos::ParameterList);
    aztecList = Teuchos::getParametersFromXmlFile("datafiles/SolverParamList.xml");
    linearSolver.setParameters(*aztecList);

    typedef LifeV::PreconditionerML         precML_type;
    typedef shp<precML_type>    precMLPtr_type;
    precML_type * precRawPtr;
    precRawPtr = new precML_type;

    GetPot dummyDatafile;
    precRawPtr->setDataFromGetPot(dummyDatafile, "precMLL");
    shp<LifeV::Preconditioner> precPtr;
    precPtr.reset(precRawPtr);

    linearSolver.setPreconditioner(precPtr);
    linearSolver.setRightHandSide(weakWSSUnique);
    linearSolver.solve(WSS);
}

void
StokesModel::
initializeVelocityFESpace(EPETRACOMM comm)
{
    // initialize fespace velocity
    M_velocityFESpace.reset(new FESPACE(this->M_treeNode->M_block->getMesh(),
                                        M_velocityOrder, 3, comm));

    M_velocityFESpaceETA.reset(new ETFESPACE3(M_velocityFESpace->mesh(),
                                              &(M_velocityFESpace->refFE()),
                                              comm));
}


void
StokesModel::
initializePressureFESpace(EPETRACOMM comm)
{
    // initialize fespace pressure
    M_pressureFESpace.reset(new FESPACE(this->M_treeNode->M_block->getMesh(),
                                        M_pressureOrder, 1, comm));
    M_pressureFESpaceETA.reset(new ETFESPACE1(M_pressureFESpace->mesh(),
                                              &(M_pressureFESpace->refFE()),
                                              comm));
}

}
