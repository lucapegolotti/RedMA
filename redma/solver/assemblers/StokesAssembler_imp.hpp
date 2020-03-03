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
    LifeV::LifeChrono chrono;
    chrono.start();

    std::string msg = "[";
    msg += this->M_name;
    msg += "] initializing ...";
    printlog(YELLOW, msg, this->M_data.getVerbose());

    initializeFEspaces();

    M_stiffness = assembleStiffness();
    M_mass = assembleMass();
    M_divergence = assembleDivergence();

    assembleFlowRateVectors();
    // assembleFlowRateJacobians();

    setExporter();

    msg = "done, in ";
    msg += std::to_string(chrono.diff());
    msg += " seconds\n";
    printlog(YELLOW, msg, this->M_data.getVerbose());
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
        std::cout << "\nInlet normal\n" << std::endl;
        face.print();
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
getNorm(const unsigned int& fieldIndex)
{
    InMatrixType retMat;
    if (fieldIndex == 0)
    {
        retMat = M_mass.block(0,0);
        retMat += M_stiffness.block(0,0);
    }
    else
        throw new Exception("Mass pressure has not been implemented yet");

    return retMat;
}

template <class InVectorType, class InMatrixType>
InMatrixType
StokesAssembler<InVectorType, InMatrixType>::
getConstraintMatrix()
{
    return M_divergence.block(0,1);
}

}
