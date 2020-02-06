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

    M_useLifting = this->M_data("bc_conditions/lifting", true);
}

template <class InVectorType, class InMatrixType>
void
StokesAssembler<InVectorType,InMatrixType>::
setup()
{
    LifeV::LifeChrono chrono;
    chrono.start();

    printlog(YELLOW, "[StokesAssembler] initializing "
                     "StokesAssembler ...", this->M_data.getVerbose());

    initializeFEspaces();

    assembleStiffness();
    assembleMass();
    assembleDivergence();

    assembleFlowRateVectors();
    assembleFlowRateJacobians();

    setExporter();

    if (M_useLifting)
    {
        M_TMAlifting = TimeMarchingAlgorithmFactory<InVectorType COMMA InMatrixType>(this->M_data);
        M_TMAlifting->setup(getZeroVector());
    }

    std::string msg = "done, in ";
    msg += std::to_string(chrono.diff());
    msg += " seconds\n";
    printlog(YELLOW, msg, this->M_data.getVerbose());
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
exportSolution(const double& t, const BlockVector<InVectorType>& sol)
{
    if (M_useLifting)
    {
        BlockVector<InVectorType> copySol;
        copySol.hardCopy(sol);

        this->M_bcManager->applyDirichletBCs(t, copySol, getFESpaceBCs(),
                                             getComponentBCs());
        *M_velocityExporter = *copySol.block(0).data();
        *M_pressureExporter = *copySol.block(1).data();
    }
    else
    {
        *M_velocityExporter = *sol.block(0).data();
        *M_pressureExporter = *sol.block(1).data();
    }

    BlockVector<InVectorType> solCopy(2);
    solCopy.block(0).data() = M_velocityExporter;
    computeFlowRates(solCopy, true);

    CoutRedirecter ct;
    ct.redirect();
    M_exporter->postProcess(t);
    printlog(CYAN, ct.restore());
}

template <class InVectorType, class InMatrixType>
void
StokesAssembler<InVectorType, InMatrixType>::
postProcess(const double& t, const BlockVector<InVectorType>& sol)
{
    // shift solutions in multistep method embedded in windkessels
    this->M_bcManager->postProcess();

    // shift liftings for derivative approximation
    if (M_useLifting)
    {
        BlockVector<InVectorType> lifting = computeLifting(t);
        M_TMAlifting->shiftSolutions(lifting);
    }
}

template <class InVectorType, class InMatrixType>
BlockMatrix<InMatrixType>
StokesAssembler<InVectorType, InMatrixType>::
getMass(const double& time, const BlockVector<InVectorType>& sol)
{
    return M_mass;
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

    BlockVector<InVectorType> lifting;
    BlockVector<InVectorType> liftingDt;
    if (M_useLifting)
    {
        lifting = computeLifting(time);
        retVec += (systemMatrix * lifting);

        double dt = this->M_data("time_discretization/dt", 0.01);
        liftingDt = M_TMAlifting->computeDerivative(lifting, dt);
        retVec -= (M_mass * liftingDt);
    }

    addNeumannBCs(retVec, time, sol);

    apply0DirichletBCs(retVec);

    return retVec;
}

template <class InVectorType, class InMatrixType>
void
StokesAssembler<InVectorType, InMatrixType>::
apply0DirichletBCs(BlockVector<InVectorType>& vector)
{
    if (M_useLifting)
        this->M_bcManager->apply0DirichletBCs(vector, getFESpaceBCs(), getComponentBCs());
}


template <class InVectorType, class InMatrixType>
BlockVector<FEVECTOR>
StokesAssembler<InVectorType, InMatrixType>::
computeLifting(const double& time) const
{
    BlockVector<InVectorType> lifting;
    lifting.resize(2);
    lifting.block(0).data().reset(new VECTOREPETRA(M_velocityFESpace->map(),
                                                   LifeV::Unique));
    lifting.block(0).data()->zero();
    lifting.block(1).data().reset(new VECTOREPETRA(M_pressureFESpace->map(),
                                                   LifeV::Unique));
    lifting.block(1).data()->zero();

    this->M_bcManager->applyDirichletBCs(time, lifting, getFESpaceBCs(),
                                         getComponentBCs());
    return lifting;
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

        flowRates[face.M_flag] = std::abs(sol.block(0).data()->dot(*M_flowRateVectors[face.M_flag]));
        msg = "[StokesAssembler] inflow rate = ";
        msg += std::to_string(flowRates[face.M_flag]);
        msg += "\n";
        printlog(YELLOW, msg, verbose);
    }

    if (this->M_treeNode->isOutletNode())
    {
        auto faces = this->M_treeNode->M_block->getOutlets();

        for (auto face : faces)
        {
            flowRates[face.M_flag] = std::abs(sol.block(0).data()->dot(*M_flowRateVectors[face.M_flag]));
            msg = "[StokesAssembler] outflow rate = ";
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
    auto flowRates = computeFlowRates(sol);

    if (this->M_treeNode->isOutletNode())
        this->M_bcManager->applyNeumannBc(time, input, M_velocityFESpace, 0, flowRates);
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

    if (this->M_treeNode->isOutletNode())
    {
        auto flowRates = computeFlowRates(sol);

        for (auto rate : flowRates)
        {
            double dhdQ = this->M_bcManager->getNeumannJacobian(time, rate.first, rate.second);
            BlockMatrix<InMatrixType> curjac(this->M_nComponents, this->M_nComponents);

            curjac.block(0,0).data().reset(new MATRIXEPETRA(*M_flowRateJacobians[rate.first]));
            curjac.block(0,0) *= dhdQ;

            retMat += curjac;
        }
    }

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

    M_exporter->setPostDir(outdir);
}

}
