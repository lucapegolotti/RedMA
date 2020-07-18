namespace RedMA
{

template <class InVectorType, class InMatrixType>
StokesAssembler<InVectorType, InMatrixType>::
StokesAssembler(const DataContainer& data,
                SHP(TreeNode) treeNode) :
  aAssembler<InVectorType, InMatrixType>(data, treeNode)
{
    M_density = this->M_data("fluid/density", 1.0);
    M_viscosity = this->M_data("fluid/viscosity", 0.035);
    this->M_nComponents = 2;
    this->M_bcManager.reset(new BCManager(data, treeNode));
    this->M_name = "StokesAssembler";
}

template <class InVectorType, class InMatrixType>
void
StokesAssembler<InVectorType,InMatrixType>::
setup()
{
    Chrono chrono;
    chrono.start();

    std::string msg = "[";
    msg += this->M_name;
    msg += "] initializing ...";
    printlog(YELLOW, msg, this->M_data.getVerbose());

    initializeFEspaces();

    M_mass = assembleMass(); // #1
    M_stiffness = assembleStiffness(); // #2
    M_divergence = assembleDivergence(); // #3

    assembleFlowRateVectors();
    // assembleFlowRateJacobians();

    setExporter();

    // initializePythonStructures();

    msg = "done, in ";
    msg += std::to_string(chrono.diff());
    msg += " seconds\n";
    printlog(YELLOW, msg, this->M_data.getVerbose());
}

template <class InVectorType, class InMatrixType>
void
StokesAssembler<InVectorType, InMatrixType>::
initializePythonStructures()
{
    // setenv("PYTHONPATH",".",1);

    // Py_Initialize();
    // PyObject* pName = PyUnicode_DecodeFSDefault("test");

    // M_pModule = PyImport_Import(pName);
    // Py_DECREF(pName);

    // M_pFunc = PyObject_GetAttrString(M_pModule, "evaluate_model");
}

template <class InVectorType, class InMatrixType>
SHP(VECTOREPETRA)
StokesAssembler<InVectorType, InMatrixType>::
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

template <class InVectorType, class InMatrixType>
BlockMatrix<InMatrixType>
StokesAssembler<InVectorType, InMatrixType>::
getMassJacobian(const double& time, const BlockVector<InVectorType>& sol)
{
    BlockMatrix<InMatrixType> retMat(this->M_nComponents, this->M_nComponents);
    return retMat;
}

template <class InVectorType, class InMatrixType>
BlockVector<InVectorType>
StokesAssembler<InVectorType, InMatrixType>::
getForcingTerm(const double& time) const
{
    // for the moment we don't consider forcing terms
    return getZeroVector();
}

template <class InVectorType, class InMatrixType>
void
StokesAssembler<InVectorType, InMatrixType>::
initializeFEspaces()
{
    // initialize fespace velocity
    std::string orderVelocity = this->M_data("fluid/velocity_order", "P2");
    M_velocityFESpace.reset(new FESPACE(this->M_treeNode->M_block->getMesh(),
                                        orderVelocity, 3, this->M_comm));

    M_velocityFESpaceETA.reset(new ETFESPACE3(M_velocityFESpace->mesh(),
                                            &(M_velocityFESpace->refFE()),
                                              this->M_comm));

    // initialize fespace velocity
    std::string orderPressure = this->M_data("fluid/pressure_order", "P1");

    M_pressureFESpace.reset(new FESPACE(this->M_treeNode->M_block->getMesh(),
                                        orderPressure, 1, this->M_comm));
    M_pressureFESpaceETA.reset(new ETFESPACE1(M_pressureFESpace->mesh(),
                                            &(M_pressureFESpace->refFE()),
                                              this->M_comm));

}

template <class InVectorType, class InMatrixType>
void
StokesAssembler<InVectorType, InMatrixType>::
postProcess(const double& t, const BlockVector<InVectorType>& sol)
{
    // shift solutions in multistep method embedded in windkessels
    this->M_bcManager->postProcess();
}

template <class InVectorType, class InMatrixType>
BlockMatrix<InMatrixType>
StokesAssembler<InVectorType, InMatrixType>::
getMass(const double& time, const BlockVector<InVectorType>& sol)
{
    return M_mass;
}

template <class InVectorType, class InMatrixType>
void
StokesAssembler<InVectorType, InMatrixType>::
assembleFlowRateVectors()
{
    // assemble inflow flow rate vector
    if (this->M_treeNode->isInletNode())
    {
        auto face = this->M_treeNode->M_block->getInlet();
        // std::cout << "\nInlet normal\n" << std::endl;
        // face.print();
        this->M_flowRateVectors[face.M_flag] = assembleFlowRateVector(face);
    }

    if (this->M_treeNode->isOutletNode())
    {
        auto faces = this->M_treeNode->M_block->getOutlets();

        for (auto face : faces)
            this->M_flowRateVectors[face.M_flag] = assembleFlowRateVector(face);
    }
}

template <class InVectorType, class InMatrixType>
BlockVector<InVectorType>
StokesAssembler<InVectorType, InMatrixType>::
getRightHandSide(const double& time, const BlockVector<InVectorType>& sol)
{
    BlockVector<InVectorType> retVec;
    BlockMatrix<InMatrixType> systemMatrix;

    systemMatrix.resize(this->M_nComponents, this->M_nComponents);
    systemMatrix += M_stiffness;
    systemMatrix += M_divergence;
    systemMatrix *= (-1.0);

    retVec.softCopy(systemMatrix * sol);

    // addNeumannBCs(retVec, time, sol);

    return retVec;
}

template <class InVectorType, class InMatrixType>
std::map<unsigned int, double>
StokesAssembler<InVectorType, InMatrixType>::
computeFlowRates(const BlockVector<InVectorType>& sol, bool verbose)
{
    std::string msg;
    std::map<unsigned int, double> flowRates;
    if (this->M_treeNode->isInletNode())
    {
        auto face = this->M_treeNode->M_block->getInlet();

        flowRates[face.M_flag] = sol.block(0).data()->dot(*M_flowRateVectors[face.M_flag]);
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
            flowRates[face.M_flag] = sol.block(0).data()->dot(*M_flowRateVectors[face.M_flag]);
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


template <class InVectorType, class InMatrixType>
void
StokesAssembler<InVectorType, InMatrixType>::
addNeumannBCs(BlockVector<FEVECTOR>& input, const double& time,
              const BlockVector<InVectorType>& sol)
{
    // if (this->M_treeNode->isOutletNode())
    //     this->M_bcManager->applyNeumannBc(time, input, M_velocityFESpace, 0, flowRates);

    if (this->M_treeNode->isOutletNode())
    {
        auto flowRates = computeFlowRates(sol, true);

        for (auto rate : flowRates)
        {
            double P = this->M_bcManager->getNeumannBc(time, rate.first, rate.second);
            BlockVector<InVectorType> curvec(this->M_nComponents);

            curvec.block(0).data().reset(new VECTOREPETRA(*M_flowRateVectors[rate.first]));
            curvec.block(0) *= P;
            input += curvec;

            // addBackFlowStabilization(input, sol, rate.first);
        }
    }
}

template <class InVectorType, class InMatrixType>
void
StokesAssembler<InVectorType, InMatrixType>::
computeWallShearStress(SHP(VECTOREPETRA) velocity)
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

    M_WSSExporter->zero();

    SHP(VECTOREPETRA) velocityRepeated(new VECTOREPETRA(*velocity,
                                                         Repeated));
    SHP(VECTOREPETRA) weakWSS(new VECTOREPETRA(M_velocityFESpace->map()));

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
             ) >> weakWSS;


    M_WSSExporter->zero();

    LinearSolver linearSolver(this->M_comm);
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
    linearSolver.setRightHandSide(weakWSS);
    linearSolver.solve(M_WSSExporter);
}

template <class InVectorType, class InMatrixType>
BlockMatrix<InMatrixType>
StokesAssembler<InVectorType, InMatrixType>::
getJacobianRightHandSide(const double& time, const BlockVector<InVectorType>& sol)
{
    BlockMatrix<InMatrixType> retMat;
    retMat.resize(this->M_nComponents, this->M_nComponents);

    retMat += M_stiffness;
    retMat += M_divergence;

    retMat *= (-1.0);

    // ATTENTION: here I should add the part relative to Neumann conditions
    // if they depend on the solution (as with 0D coupling)

    // if (this->M_treeNode->isOutletNode())
    // {
    //     auto flowRates = computeFlowRates(sol);
    //
    //     for (auto rate : flowRates)
    //     {
    //         double dhdQ = this->M_bcManager->getNeumannJacobian(time, rate.first, rate.second);
    //         BlockMatrix<InMatrixType> curjac;
    //         curjac.hardCopy(M_flowRateJacobians[rate.first]);
    //         curjac *= dhdQ;
    //
    //         retMat += curjac;
    //     }
    // }

    // this->M_bcManager->apply0DirichletMatrix(retMat, getFESpaceBCs(),
    //                                    getComponentBCs(), 0.0);

    return retMat;
}

template <class InVectorType, class InMatrixType>
void
StokesAssembler<InVectorType, InMatrixType>::
setExporter()
{
    std::string outputName = "block";
    outputName += std::to_string(this->M_treeNode->M_ID);

    std::string outdir = this->M_data("exporter/outdir", "solutions/");
    boost::filesystem::create_directory(outdir);

    std::string format = this->M_data("exporter/type", "hdf5");
    if (!std::strcmp(format.c_str(), "hdf5"))
        M_exporter.reset(new LifeV::ExporterHDF5<MESH>(this->M_data.getDatafile(),
                                                       outputName));
    else
        M_exporter.reset(new LifeV::ExporterVTK<MESH>(this->M_data.getDatafile(),
                                                      outputName));
    M_exporter->setMeshProcId(M_velocityFESpace->mesh(), this->M_comm->MyPID());

    M_velocityExporter.reset(new VECTOREPETRA(M_velocityFESpace->map(),
                                              M_exporter->mapType()));

    M_pressureExporter.reset(new VECTOREPETRA(M_pressureFESpace->map(),
                                              M_exporter->mapType()));

    M_exporter->addVariable(LifeV::ExporterData<MESH>::VectorField,
                         "velocity", M_velocityFESpace, M_velocityExporter, 0.0);

    M_exporter->addVariable(LifeV::ExporterData<MESH>::ScalarField,
                         "pressure", M_pressureFESpace, M_pressureExporter, 0.0);

    bool exportWSS = this->M_data("exporter/export_wss", 1);
    if (exportWSS)
    {
        M_WSSExporter.reset(new VECTOREPETRA(M_velocityFESpace->map(),
                                             M_exporter->mapType()));

        M_exporter->addVariable(LifeV::ExporterData<MESH>::VectorField,
                                "WSS", M_velocityFESpace, M_WSSExporter, 0.0);
    }

    M_exporter->setPostDir(outdir);
}

template <class InVectorType, class InMatrixType>
SHP(FESPACE)
StokesAssembler<InVectorType, InMatrixType>::
getFEspace(unsigned int index) const
{
    if (index == 0)
        return M_velocityFESpace;
    else if (index == 1)
        return M_pressureFESpace;

    return nullptr;
}

template <class InVectorType, class InMatrixType>
BlockVector<FEVECTOR>
StokesAssembler<InVectorType, InMatrixType>::
getFELifting(const double& time) const
{
    BlockVector<FEVECTOR> lifting;
    lifting.resize(2);
    lifting.block(0).data().reset(new VECTOREPETRA(M_velocityFESpace->map(),
                                                   LifeV::Unique));
    lifting.block(0).data()->zero();
    lifting.block(1).data().reset(new VECTOREPETRA(M_pressureFESpace->map(),
                                                   LifeV::Unique));
    lifting.block(1).data()->zero();

    this->M_bcManager->applyDirichletBCs(time, lifting, this->getFESpaceBCs(),
                                         this->getComponentBCs());

    return lifting;
}

template <class InVectorType, class InMatrixType>
std::vector<BlockMatrix<InMatrixType>>
StokesAssembler<InVectorType, InMatrixType>::
getMatrices() const
{
    std::vector<BlockMatrix<InMatrixType>> retVec;

    retVec.push_back(M_mass);
    retVec.push_back(M_stiffness);
    retVec.push_back(M_divergence);

    return retVec;
}

template <class InVectorType, class InMatrixType>
BlockMatrix<InMatrixType>
StokesAssembler<InVectorType, InMatrixType>::
assembleMatrix(const unsigned int& index, BlockMDEIMStructure* structure)
{
    if (index == 0)
    {
        return assembleMass(structure);
    }
    else if (index == 1)
    {
        return assembleStiffness(structure);
    }
    else if (index == 2)
    {
        return assembleDivergence(structure);
    }
}

template <class InVectorType, class InMatrixType>
InMatrixType
StokesAssembler<InVectorType, InMatrixType>::
getConstraintMatrix()
{
    return M_divergence.block(0,1);
}

template <class InVectorType, class InMatrixType>
BlockMatrix<MatrixEp>
StokesAssembler<InVectorType,InMatrixType>::
assembleReducedStiffness(BlockMDEIMStructure* structure)
{
    using namespace LifeV;
    using namespace ExpressionAssembly;

    BlockMatrix<MatrixEp> stiffness;

    stiffness.resize(this->M_nComponents, this->M_nComponents);
    bool useFullStrain = this->M_data("fluid/use_strain", true);

    SHP(MatrixEpetra<double>) A(new MatrixEpetra<double>(M_velocityFESpace->map()));

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

    stiffness.block(0,0).data() = A;

    this->M_bcManager->apply0DirichletMatrix(stiffness, getFESpaceBCs(),
                                             getComponentBCs(), 0.0);

    return stiffness;
}

template <class InVectorType, class InMatrixType>
BlockMatrix<MatrixEp>
StokesAssembler<InVectorType,InMatrixType>::
assembleReducedMass(BlockMDEIMStructure* structure)
{
    using namespace LifeV;
    using namespace ExpressionAssembly;

    BlockMatrix<MatrixEp> mass;

    mass.resize(this->M_nComponents, this->M_nComponents);

    SHP(MatrixEpetra<double>) M(new MatrixEpetra<double>(M_velocityFESpace->map()));

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
    mass.block(0,0).data() = M;

    this->M_bcManager->apply0DirichletMatrix(mass, getFESpaceBCs(),
                                             getComponentBCs(), 1.0);

    return mass;
}

template <class InVectorType, class InMatrixType>
BlockMatrix<MatrixEp>
StokesAssembler<InVectorType,InMatrixType>::
assembleReducedDivergence(BlockMDEIMStructure* structure)
{
    using namespace LifeV;
    using namespace ExpressionAssembly;

    BlockMatrix<MatrixEp> divergence;

    divergence.resize(this->M_nComponents, this->M_nComponents);

    SHP(MatrixEpetra<double>) BT(new MatrixEpetra<double>(this->M_velocityFESpace->map()));

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

    SHP(MatrixEpetra<double>) B(new MatrixEpetra<double>(M_pressureFESpace->map()));

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

    divergence.block(0,1).data() = BT;
    divergence.block(1,0).data() = B;

    this->M_bcManager->apply0DirichletMatrix(divergence, getFESpaceBCs(),
                                             getComponentBCs(), 0.0);

    return divergence;
}

template <class InVectorType, class InMatrixType>
void
StokesAssembler<InVectorType,InMatrixType>::
setMDEIMs(SHP(MDEIMManager) mdeimManager)
{
    std::string mdeimdir = this->M_data("rb/online/mdeim/directory", "mdeims");
    std::string meshName = this->M_treeNode->M_block->getMeshName();
    unsigned int dashpos = meshName.find("/");
    unsigned int formatpos = meshName.find(".mesh");
    std::string actualName = meshName.substr(dashpos + 1,
                                             formatpos - dashpos - 1);

    auto mdeims = mdeimManager->getMDEIMS(actualName);

    M_mdeimMass = mdeims[0];
    M_mdeimStiffness = mdeims[1];
    M_mdeimDivergence = mdeims[2];
}

template <class InVectorType, class InMatrixType>
void
StokesAssembler<InVectorType,InMatrixType>::
setRBBases(SHP(RBBasesManager) rbManager)
{
    std::string mdeimdir = this->M_data("rb/online/mdeim/directory", "mdeims");
    std::string meshName = this->M_treeNode->M_block->getMeshName();
    unsigned int dashpos = meshName.find("/");
    unsigned int formatpos = meshName.find(".mesh");
    std::string actualName = meshName.substr(dashpos + 1,
                                             formatpos - dashpos - 1);

    // beware that at this point the rb bases have not been loaded yet
    M_bases = rbManager->getRBBases(actualName);
    M_bases->setFESpace(M_velocityFESpace, 0);
    M_bases->setFESpace(M_pressureFESpace, 1);
}

template <class InVectorType, class InMatrixType>
BlockVector<FEVECTOR>
StokesAssembler<InVectorType,InMatrixType>::
convertFunctionRBtoFEM(BlockVector<RBVECTOR> rbSolution) const
{
    BlockVector<FEVECTOR> retVec(2);

    unsigned int id = this->M_treeNode->M_ID;

    if (rbSolution.block(0).data())
        retVec.block(0).data() = M_bases->reconstructFEFunction(rbSolution.block(0), 0, id);

    if (rbSolution.block(1).data())
        retVec.block(1).data() = M_bases->reconstructFEFunction(rbSolution.block(1), 1, id);

    return retVec;
}

template <class InVectorType, class InMatrixType>
void
StokesAssembler<InVectorType,InMatrixType>::
exportNorms(double t)
{
    bool exportNorms = this->M_data("exporter/exportnorms", true);

    if (exportNorms)
    {
        std::string outputName = this->M_data("exporter/outdir", "solutions/") + "block";
        outputName += std::to_string(this->M_treeNode->M_ID);
        outputName += "_norm.txt";
        std::ofstream filename(outputName, std::ios_base::app);
        filename << t << ",";
        filename << M_velocityFESpace->h1Norm(*M_velocityExporter) << ",";
        filename << M_pressureFESpace->l2Norm(*M_pressureExporter) << "\n";
    }
}

template <class InVectorType, class InMatrixType>
void
StokesAssembler<InVectorType,InMatrixType>::
setExtrapolatedSolution(const BlockVector<InVectorType>& exSol)
{
    M_extrapolatedSolution.softCopy(exSol);
}

template <class InVectorType, class InMatrixType>
void
StokesAssembler<InVectorType,InMatrixType>::
applyPiola(BlockVector<FEVECTOR> solution, bool inverse)
{
    std::string msg = "[";
    msg += this->M_name;
    msg += "] apply Piola\n";
    printlog(YELLOW, msg, this->M_data.getVerbose());

    using namespace LifeV;
    using namespace ExpressionAssembly;

    auto defAssembler = this->M_defaultAssemblers->
                        getDefaultAssembler(this->M_treeNode->M_block->getMeshName());

    // auto refDivergence = defAssembler->assembleMatrix(2);
    // SHP(MatrixEpetra<double>) B(new MatrixEpetra<double>(M_pressureFESpace->map()));
    //
    // integrate(elements(M_velocityFESpaceETA->mesh()),
    //          M_pressureFESpace->qr(),
    //          M_pressureFESpaceETA,
    //          M_velocityFESpaceETA,
    //          phi_i * div(phi_j)
    //      ) >> B;
    //
    // B->globalAssemble(M_velocityFESpace->mapPtr(),
    //                   M_pressureFESpace->mapPtr());
    //
    // FEMATRIX Bwrap;
    // Bwrap.data() = B;
    //
    // FEVECTOR res1 = Bwrap * solution.block(0);
    //
    // std::cout << this->M_treeNode->M_block->getMeshName() << std::endl << std::flush;
    // std::cout << "norm 1 = " << res1.norm2() << std::endl << std::flush;
    //
    // FEVECTOR res2 = refDivergence.block(1,0) * solution.block(0);
    //
    // std::cout << "norm 2 = " << res2.norm2() << std::endl << std::flush;

    auto defVelocityFESpace = defAssembler->getFEspace(0);

    auto velocity = solution.block(0).data();

    Epetra_Map epetraMap = velocity->epetraMap();
    unsigned int numElements = epetraMap.NumMyElements() / 3;

    // find xs ys and zs of mesh
    if (M_xs == nullptr)
    {
        M_xs.reset(new VECTOREPETRA(defVelocityFESpace->map()));

        defVelocityFESpace->interpolate([](const double& t,
                                           const double& x,
                                           const double& y,
                                           const double& z,
                                           const unsigned int & i) {return x;},
                                           *M_xs, 0.0);
    }

    if (M_ys == nullptr)
    {
        M_ys.reset(new VECTOREPETRA(defVelocityFESpace->map()));

        defVelocityFESpace->interpolate([](const double& t,
                                           const double& x,
                                           const double& y,
                                           const double& z,
                                           const unsigned int & i) {return y;},
                                           *M_ys, 0.0);
    }

    if (M_zs == nullptr)
    {
        M_zs.reset(new VECTOREPETRA(defVelocityFESpace->map()));

        defVelocityFESpace->interpolate([](const double& t,
                                           const double& x,
                                           const double& y,
                                           const double& z,
                                           const unsigned int & i) {return z;},
                                           *M_zs, 0.0);
    }

    LifeV::MatrixSmall<3,3>* transformationMatrix;
    double determinant;
    if (inverse)
        transformationMatrix = new LifeV::MatrixSmall<3,3>();


    for (unsigned int dof = 0; dof < numElements; dof++)
    {
        double& x = M_xs->operator[](dof);
        double& y = M_ys->operator[](dof);
        double& z = M_zs->operator[](dof);

        auto jacobian = this->M_treeNode->M_block->computeJacobianGlobalTransformation(x, y, z);

        double xx = x, yy = y, zz = z;

        this->M_treeNode->M_block->globalTransf(xx, yy, zz);

        double& a = jacobian(0,0);
        double& b = jacobian(0,1);
        double& c = jacobian(0,2);
        double& d = jacobian(1,0);
        double& e = jacobian(1,1);
        double& f = jacobian(1,2);
        double& g = jacobian(2,0);
        double& h = jacobian(2,1);
        double& i = jacobian(2,2);

        determinant = std::abs(a*e*i - a*f*h - b*d*i + b*f*g + c*d*h - c*e*g);

        if (inverse)
        {
            this->M_treeNode->M_block->matrixInverse(jacobian, transformationMatrix);
            determinant = 1.0/determinant;
        }
        else
            transformationMatrix = &jacobian;

        LifeV::VectorSmall<3> curU;
        curU(0) = velocity->operator[](dof);
        curU(1) = velocity->operator[](dof + numElements);
        curU(2) = velocity->operator[](dof + numElements * 2);

        LifeV::VectorSmall<3> res;
        res = 1./determinant * (*transformationMatrix) * curU;

        velocity->operator[](dof) = res(0);
        velocity->operator[](dof + numElements) = res(1);
        velocity->operator[](dof + numElements * 2) = res(2);
    }

    // res1 = Bwrap * solution.block(0);
    // std::cout << "norm 1 = " << res1.norm2() << std::endl << std::flush;
    // res2 = refDivergence.block(1,0) * solution.block(0);
    // std::cout << "norm 2 = " << res2.norm2() << std::endl << std::flush;

    if (inverse)
        delete transformationMatrix;

    // defAssembler->exportSolution(0.0, solution);
}

}
