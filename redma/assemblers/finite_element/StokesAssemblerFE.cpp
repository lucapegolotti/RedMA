#include "StokesAssemblerFE.hpp"

namespace RedMA
{

StokesAssemblerFE::
StokesAssemblerFE(const DataContainer& data,
                  shp<TreeNode> treeNode) :
  aAssemblerFE(data, treeNode)
{
    this->M_bcManager.reset(new BCManager(data, treeNode));
    M_name = "StokesAssemblerFE";
    M_nComponents = 2;
    M_density = data("fluid/density", 1.0);
    M_viscosity = data("fluid/viscosity", 0.035);
    M_velocityOrder = data("fluid/velocity_order", "P2");
    M_pressureOrder = data("fluid/pressure_order", "P1");
}

void
StokesAssemblerFE::
setup()
{
    Chrono chrono;
    chrono.start();

    std::string msg = "[";
    msg += this->M_name;
    msg += "] initializing ...";
    printlog(YELLOW, msg, this->M_data.getVerbose());
    initializeFEspaces();
    M_mass = spcast<BlockMatrix>(assembleMass(M_bcManager)); // #1
    M_stiffness = spcast<BlockMatrix>(assembleStiffness(M_bcManager)); // #2
    M_divergence = spcast<BlockMatrix>(assembleDivergence(M_bcManager)); // #3
    assembleFlowRateVectors();
    // assembleFlowRateJacobians(this->M_bcManager);

    setExporter();

    // initializePythonStructures();

    msg = "done, in ";
    msg += std::to_string(chrono.diff());
    msg += " seconds\n";
    printlog(YELLOW, msg, this->M_data.getVerbose());
}

void
StokesAssemblerFE::
postProcess(const double& t,
            const shp<aVector>& sol)
{
    // shift solutions in multistep method embedded in windkessels
    M_bcManager->postProcess();
}


std::map<unsigned int, std::vector<shp<BlockVector>>>
StokesAssemblerFE::
importSolution(const std::string& filename) const
{
    std::fstream inVel(filename + "velocity.txt");
    std::fstream inPres(filename + "pressure.txt");
    std::string line;
    std::vector<shp<VECTOREPETRA>> vecVelocity;
    std::vector<shp<VECTOREPETRA>> vecPressure;

    std::map<unsigned int, std::vector<shp<BlockVector>>> retMap;

    unsigned int cnt = 0;
    std::vector<LifeV::Int> indicesPres;

    while(std::getline(inPres, line))
    {
        shp<VECTOREPETRA> tmpEpetraVecPressure (new VECTOREPETRA(M_pressureFESpace->map(),
                                                                 LifeV::Unique));

        double value;
        std::stringstream ss(line);

        std::vector<double> values;
        while (ss >> value)
            values.push_back(value);

        if (cnt == 0)
            for (LifeV::Int i=0; i<tmpEpetraVecPressure->size(); ++i)
                indicesPres.push_back(i);
        cnt += 1;

        tmpEpetraVecPressure->setCoefficients(indicesPres, values);

        vecPressure.push_back(tmpEpetraVecPressure);
    }


    cnt = 0;
    std::vector<LifeV::Int> indicesVel;

    while (std::getline(inVel, line))
    {
        shp<VECTOREPETRA> tmpEpetraVecVelocity (new VECTOREPETRA(M_velocityFESpace->map(),
                                                                 LifeV::Unique));

        double value;
        std::stringstream ss(line);

        std::vector<double> values;
        while (ss >> value)
            values.push_back(value);

        if (cnt == 0)
            for (LifeV::Int i=0; i<tmpEpetraVecVelocity->size(); ++i)
                indicesVel.push_back(i);
        cnt+=1;

        tmpEpetraVecVelocity->setCoefficients(indicesVel, values);

        vecVelocity.push_back(tmpEpetraVecVelocity);
    }

    if ((vecVelocity.size() == 0) || (vecPressure.size() == 0))
        throw new Exception("Importing error! Impossible to load velocity and pressure");

    for (unsigned int count = 0; count<vecVelocity.size(); ++count)
    {
        shp<BlockVector> retBlock (new BlockVector(2));
        retBlock->setBlock(0, wrap(vecVelocity[count]));
        retBlock->setBlock(1, wrap(vecPressure[count]));
        retMap[0].push_back(retBlock);
    }

    return retMap;
}

void
StokesAssemblerFE::
exportSolution(const double& t,
               const shp<aVector>& sol)
{
    auto solBlck = convert<BlockVector>(sol);
    *M_velocityExporter = *spcast<VECTOREPETRA>(solBlck->block(0)->data());
    *M_pressureExporter = *spcast<VECTOREPETRA>(solBlck->block(1)->data());

    bool exportWSS = this->M_data("exporter/export_wss", 1);
    if (exportWSS)
        computeWallShearStress(M_velocityExporter, M_WSSExporter, M_comm);

    shp<BlockVector> solCopy(new BlockVector(2));
    solCopy->block(0)->setData(M_velocityExporter);
    computeFlowRates(solCopy, true);

    exportNorms(t, M_velocityExporter, M_pressureExporter);

    CoutRedirecter ct;
    ct.redirect();
    M_exporter->postProcess(t);
    printlog(CYAN, ct.restore());
}

shp<aMatrix>
StokesAssemblerFE::
getMass(const double& time,
        const shp<aVector>& sol)
{
    return M_mass;
}

shp<aMatrix>
StokesAssemblerFE::
getMassJacobian(const double& time,
                const shp<aVector>& sol)
{
    shp<aMatrix> retMat(new BlockMatrix(this->M_nComponents,
                                        this->M_nComponents));
    return retMat;
}

shp<aVector>
StokesAssemblerFE::
getRightHandSide(const double& time,
                 const shp<aVector>& sol)
{
    shp<BlockMatrix> systemMatrix(new BlockMatrix(this->M_nComponents,
                                                  this->M_nComponents));
    systemMatrix->add(M_stiffness);
    systemMatrix->add(M_divergence);
    systemMatrix->multiplyByScalar(-1.0);

    shp<aVector> retVec = systemMatrix->multiplyByVector(sol);
    addNeumannBCs(time, sol, retVec);

    this->M_bcManager->apply0DirichletBCs(*spcast<BlockVector>(retVec),
                                          this->getFESpaceBCs(),
                                          this->getComponentBCs());

    return retVec;
}

void
StokesAssemblerFE::
addNeumannBCs(double time,
              shp<aVector> sol,
              shp<aVector> rhs)
{

    if (aAssembler::M_treeNode->isOutletNode())
    {
        auto flowRates = computeFlowRates(sol, true);

        for (auto rate : flowRates)
        {
            double P = this->M_bcManager->getNeumannBc(time,
                                                       rate.first,
                                                       rate.second);

            shp<VECTOREPETRA> flowRateCopy(new VECTOREPETRA(*M_flowRateVectors[rate.first]));
            *flowRateCopy *= P;

            *spcast<VECTOREPETRA>(convert<BlockVector>(rhs)->block(0)->data()) += *flowRateCopy;
            // addBackFlowStabilization(input, sol, rate.first);
        }
    }
}

shp<aMatrix>
StokesAssemblerFE::
getJacobianRightHandSide(const double& time,
                         const shp<aVector>& sol)
{
    shp<BlockMatrix> retMat(new BlockMatrix(this->M_nComponents,
                                            this->M_nComponents));

    retMat->add(M_stiffness);
    retMat->add(M_divergence);
    retMat->multiplyByScalar(-1.0);


    if (aAssembler::M_treeNode->isOutletNode())
    {
        auto flowRates = computeFlowRates(sol);

        for (auto rate : flowRates)
        {
            double dhdQ = this->M_bcManager->getNeumannJacobian(time,
                                                                rate.first,
                                                                rate.second);
            shp<BlockMatrix> curjac(new BlockMatrix(2,2));
            curjac->deepCopy(M_flowRateJacobians[rate.first]);
            curjac->multiplyByScalar(dhdQ);

            retMat->add(curjac);
        }
    }

    this->M_bcManager->apply0DirichletMatrix(*retMat,
                                             getFESpaceBCs(),
                                             getComponentBCs(),
                                             0.0);

    return retMat;
}

shp<aVector>
StokesAssemblerFE::
getZeroVector() const
{
    return buildZeroVector();
}

shp<aVector>
StokesAssemblerFE::
getFELifting(const double& time) const
{
    shp<BlockVector> lifting(new BlockVector(2));
    // velocity
    shp<DistributedVector> ucomp(new DistributedVector());
    ucomp->setData(shp<VECTOREPETRA>(new VECTOREPETRA(M_velocityFESpace->map(),
                                                      LifeV::Unique)));
    ucomp->multiplyByScalar(0);

    shp<DistributedVector> pcomp(new DistributedVector());
    pcomp->setData(shp<VECTOREPETRA>(new VECTOREPETRA(M_pressureFESpace->map(),
                                                      LifeV::Unique)));
    pcomp->multiplyByScalar(0);

    lifting->setBlock(0, ucomp);
    lifting->setBlock(1, pcomp);

    this->M_bcManager->applyDirichletBCs(time,
                                         *lifting,
                                         this->getFESpaceBCs(),
                                         this->getComponentBCs());

    return lifting;
}

shp<aVector>
StokesAssemblerFE::
getLifting(const double& time) const
{
    return getFELifting(time);
}

void
StokesAssemblerFE::
initializeFEspaces()
{
    initializeVelocityFESpace(M_comm);
    initializePressureFESpace(M_comm);
}

void
StokesAssemblerFE::
setExporter()
{
    std::string outputName = "block";
    outputName += std::to_string(aAssembler::M_treeNode->M_ID);

    std::string outdir = this->M_data("exporter/outdir", "solutions/");
    fs::create_directory(outdir);

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

void
StokesAssemblerFE::
applyDirichletBCsMatrix(shp<aMatrix> matrix,
                        double diagCoeff) const
{
    auto matrixConverted = spcast<BlockMatrix>(matrix);
    this->M_bcManager->apply0DirichletMatrix(*matrixConverted, getFESpaceBCs(),
                                             getComponentBCs(), diagCoeff);
}

void
StokesAssemblerFE::
apply0DirichletBCs(shp<aVector> vector) const
{
    auto vectorConverted = spcast<BlockVector>(vector);
    this->M_bcManager->apply0DirichletBCs(*vectorConverted, getFESpaceBCs(),
                                          getComponentBCs());
}

void
StokesAssemblerFE::
applyDirichletBCs(const double& time,
                  shp<aVector> vector) const
{
    auto vectorConverted = spcast<BlockVector>(vector);
    this->M_bcManager->applyDirichletBCs(time,
                                         *vectorConverted,
                                         getFESpaceBCs(),
                                         getComponentBCs());
}

shp<FESPACE>
StokesAssemblerFE::
getFEspace(unsigned int index) const
{
    if (index == 0)
        return M_velocityFESpace;
    else if (index == 1)
        return M_pressureFESpace;

    return nullptr;
}

std::vector<shp<aMatrix>>
StokesAssemblerFE::
getMatrices() const
{
    std::vector<shp<aMatrix>> retVec;

    retVec.push_back(M_mass);
    retVec.push_back(M_stiffness);
    retVec.push_back(M_divergence);

    return retVec;
}

shp<aMatrix>
StokesAssemblerFE::
assembleMatrix(const unsigned int& index)
{
    if (index == 0)
    {
        return assembleMass(M_bcManager);
    }
    else if (index == 1)
    {
        return assembleStiffness(M_bcManager);
    }
    else if (index == 2)
    {
        return assembleDivergence(M_bcManager);
    }
}

shp<aMatrix>
StokesAssemblerFE::
getNorm(const unsigned int& fieldIndex,
        bool bcs)
{
    using namespace LifeV;
    using namespace ExpressionAssembly;

    shp<SparseMatrix> retMat(new SparseMatrix);
    if (fieldIndex == 0)
    {
        // if (!M_massVelocity.data())
        // {
            shp<MATRIXEPETRA> Nu(new MATRIXEPETRA(M_velocityFESpace->map()));

            integrate(elements(M_velocityFESpaceETA->mesh()),
                      M_velocityFESpace->qr(),
                      M_velocityFESpaceETA,
                      M_velocityFESpaceETA,
                      dot(phi_i, phi_j) +
                      dot(grad(phi_i),grad(phi_j))
                  ) >> Nu;

            Nu->globalAssemble();

            if (bcs)
            {
                shp<BlockMatrix> normWrap(new BlockMatrix(1,1));
                normWrap->setBlock(0,0,wrap(Nu));

                // note. Applying bcs does not change the norm if Dirichlet bcs are
                // homogeneous (=> lifting) or imposed weakly. Here we impose bcs
                // in order to have the correct conditions in the computation of the
                // supremizers (we have to solve a linear system..)
                M_bcManager->apply0DirichletMatrix(*normWrap, M_velocityFESpace, 0, 1.0);
            }

            retMat->setMatrix(Nu);
    }
    else
    {
            shp<MATRIXEPETRA> Np(new MATRIXEPETRA(M_pressureFESpace->map()));

            integrate(elements(M_pressureFESpaceETA->mesh()),
                      M_pressureFESpace->qr(),
                      M_pressureFESpaceETA,
                      M_pressureFESpaceETA,
                      phi_i * phi_j
                  ) >> Np;

            Np->globalAssemble();

            retMat->setMatrix(Np);
    }

    return retMat;
}

shp<aMatrix>
StokesAssemblerFE::
getConstraintMatrix()
{
    return M_divergence->block(0,1);
}

void
StokesAssemblerFE::
setExtrapolatedSolution(const shp<aVector>& exSol)
{
    M_extrapolatedSolution->shallowCopy(exSol);
}

void
StokesAssemblerFE::
applyPiola(shp<aVector> solution,
           bool inverse)
{
    using namespace LifeV;
    using namespace ExpressionAssembly;

    auto defAssembler = this->M_defaultAssemblers->
                        getDefaultAssembler(M_treeNode->M_block->getMeshName());

    if (defAssembler)
    {
        auto defVelocityFESpace = defAssembler->getFEspace(0);

        auto velocity = spcast<VECTOREPETRA>(convert<BlockVector>(solution)->block(0)->data());

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

            auto jacobian = aAssembler::M_treeNode->M_block->computeJacobianGlobalTransformation(x, y, z);

            double xx = x, yy = y, zz = z;

            aAssembler::M_treeNode->M_block->globalTransf(xx, yy, zz);

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
                aAssembler::M_treeNode->M_block->matrixInverse(jacobian,
                                                               transformationMatrix);
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

        if (inverse)
            delete transformationMatrix;
    }
}

shp<aVector>
StokesAssemblerFE::
getForcingTerm(const double& time) const
{
    // we don't consider forcing terms for the moment
    return buildZeroVector();
}


shp<aMatrix>
StokesAssemblerFE::
assembleStiffness(shp<BCManager> bcManager)
{
    using namespace LifeV;
    using namespace ExpressionAssembly;

    shp<BlockMatrix> stiffness(new BlockMatrix(2,2));

    bool useFullStrain = M_data("fluid/use_strain", true);

    shp<MATRIXEPETRA> A(new MATRIXEPETRA(M_velocityFESpace->map()));

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
    A->globalAssemble();

    shp<SparseMatrix> Awrapper(new SparseMatrix);
    Awrapper->setData(A);
    stiffness->setBlock(0,0,Awrapper);

    bcManager->apply0DirichletMatrix(*stiffness, M_velocityFESpace, 0, 0.0);
    return stiffness;
}

shp<aMatrix>
StokesAssemblerFE::
assembleMass(shp<BCManager> bcManager)
{
    using namespace LifeV;
    using namespace ExpressionAssembly;

    shp<BlockMatrix> mass(new BlockMatrix(2,2));
    shp<MATRIXEPETRA> M(new MATRIXEPETRA(M_velocityFESpace->map()));

    integrate(elements(M_velocityFESpaceETA->mesh()),
          M_velocityFESpace->qr(),
          M_velocityFESpaceETA,
          M_velocityFESpaceETA,
          value(M_density) * dot(phi_i, phi_j)
    ) >> M;
    M->globalAssemble();
    shp<SparseMatrix> Mwrapper(new SparseMatrix);
    Mwrapper->setData(M);
    mass->setBlock(0,0,Mwrapper);
    bcManager->apply0DirichletMatrix(*mass, M_velocityFESpace, 0, 1.0);

    return mass;
}

shp<aMatrix>
StokesAssemblerFE::
assembleDivergence(shp<BCManager> bcManager)
{
    using namespace LifeV;
    using namespace ExpressionAssembly;

    shp<BlockMatrix> divergence(new BlockMatrix(2,2));

    shp<MATRIXEPETRA> BT(new MATRIXEPETRA(this->M_velocityFESpace->map()));

    integrate(elements(M_velocityFESpaceETA->mesh()),
              M_velocityFESpace->qr(),
              M_velocityFESpaceETA,
              M_pressureFESpaceETA,
              value(-1.0) * phi_j * div(phi_i)
          ) >> BT;

    BT->globalAssemble(M_pressureFESpace->mapPtr(),
                       M_velocityFESpace->mapPtr());

    shp<MATRIXEPETRA> B(new MATRIXEPETRA(M_pressureFESpace->map()));

    integrate(elements(M_velocityFESpaceETA->mesh()),
             M_pressureFESpace->qr(),
             M_pressureFESpaceETA,
             M_velocityFESpaceETA,
             phi_i * div(phi_j)
         ) >> B;
    B->globalAssemble(M_velocityFESpace->mapPtr(),
                      M_pressureFESpace->mapPtr());

    shp<SparseMatrix> BTwrapper(new SparseMatrix);
    BTwrapper->setData(BT);

    shp<SparseMatrix> Bwrapper(new SparseMatrix);
    Bwrapper->setData(B);

    divergence->setBlock(0,1,BTwrapper);
    divergence->setBlock(1,0,Bwrapper);

    bcManager->apply0DirichletMatrix(*divergence, M_velocityFESpace, 0, 0.0);

    return divergence;
}

void
StokesAssemblerFE::
assembleFlowRateVectors()
{
    // assemble inflow flow rate vector
    if (M_treeNode->isInletNode())
    {
        auto faces = M_treeNode->M_block->getInlets();

        for (auto face : faces)
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
StokesAssemblerFE::
assembleFlowRateJacobians(shp<BCManager> bcManager)
{
    if (M_treeNode->isOutletNode())
    {
        auto faces = M_treeNode->M_block->getOutlets();

        for (auto face : faces)
        {
            shp<BlockMatrix> newJacobian(new BlockMatrix(2,2));
            shp<SparseMatrix> jacWrapper(new SparseMatrix());
            jacWrapper->setMatrix(assembleFlowRateJacobian(face));

            newJacobian->setBlock(0,0,jacWrapper);

            applyDirichletBCsMatrix(bcManager,newJacobian, 0.0);

            M_flowRateJacobians[face.M_flag] = newJacobian;

        }
    }
}

void
StokesAssemblerFE::
applyDirichletBCsMatrix(shp<BCManager> bcManager,
                        shp<aMatrix> matrix,
                        double diagCoeff)
{
    auto matrixConverted = spcast<BlockMatrix>(matrix);
    bcManager->apply0DirichletMatrix(*matrixConverted,
                                     M_velocityFESpace,
                                     0,
                                     diagCoeff);
}

std::map<unsigned int, double>
StokesAssemblerFE::
computeFlowRates(shp<aVector> sol, bool verbose)
{
    auto solBlck = convert<BlockVector>(sol);

    std::string msg;
    std::map<unsigned int, double> flowRates;
    if (M_treeNode->isInletNode())
    {
        auto faces = M_treeNode->M_block->getInlets();

        for (auto face : faces)
        {
            flowRates[face.M_flag] = spcast<VECTOREPETRA>(solBlck->block(0)->data())->dot(*M_flowRateVectors[face.M_flag]);
            std::string msg = "[";
            msg += "StokesAssemblerFE";
            msg += "]  inflow rate = ";
            msg += std::to_string(flowRates[face.M_flag]);
            msg += "\n";
            printlog(YELLOW, msg, verbose);
        }
    }

    if (this->M_treeNode->isOutletNode())
    {
        auto faces = this->M_treeNode->M_block->getOutlets();

        for (auto face : faces)
        {
            flowRates[face.M_flag] = spcast<VECTOREPETRA>(solBlck->block(0)->data())->dot(*M_flowRateVectors[face.M_flag]);
            std::string msg = "[";
            msg += "StokesAssemblerFE";
            msg += "]  outflow rate = ";
            msg += std::to_string(flowRates[face.M_flag]);
            msg += "\n";
            printlog(YELLOW, msg, verbose);
        }
    }

    return flowRates;
}

shp<VECTOREPETRA>
StokesAssemblerFE::
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
StokesAssemblerFE::
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

    // this should be optimized
    shp<MATRIXEPETRA> flowRateJacobian;
    flowRateJacobian.reset(new MATRIXEPETRA(M_velocityFESpace->map(),
                                            numGlobalElements,
                                            false));

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
StokesAssemblerFE::
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
StokesAssemblerFE::
exportNorms(double t, shp<VECTOREPETRA> velocity, shp<VECTOREPETRA> pressure)
{
    bool exportNorms = M_data("exporter/exportnorms", true);

    if (exportNorms)
    {
        std::string outputName = M_data("exporter/outdir", "solutions/") + "block";
        outputName += std::to_string(this->M_treeNode->M_ID);
        outputName += "_norm.txt";
        std::ofstream filename(outputName, std::ios_base::app);
        filename << t << ",";
        filename << M_velocityFESpace->h1Norm(*velocity) << ",";
        filename << M_pressureFESpace->l2Norm(*pressure) << "\n";
    }
}

// void
// StokesAssemblerFE::
// initializePythonStructures()
// {
//     // setenv("PYTHONPATH",".",1);
//
//     // Py_Initialize();
//     // PyObject* pName = PyUnicode_DecodeFSDefault("test");
//
//     // M_pModule = PyImport_Import(pName);
//     // Py_DECREF(pName);
//
//     // M_pFunc = PyObject_GetAttrString(M_pModule, "evaluate_model");
// }

void
StokesAssemblerFE::
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
    shp<VECTOREPETRA> weakWSSRepeated(new VECTOREPETRA(M_velocityFESpace->map(),
                                                       Repeated));

    integrate(boundary(M_velocityFESpaceETA->mesh(), wallFlag),
              myBDQR,
              M_velocityFESpaceETA,
              value(M_viscosity) *
              dot(
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
StokesAssemblerFE::
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
StokesAssemblerFE::
initializePressureFESpace(EPETRACOMM comm)
{
    // initialize fespace pressure
    M_pressureFESpace.reset(new FESPACE(this->M_treeNode->M_block->getMesh(),
                                       M_pressureOrder, 1, comm));
    M_pressureFESpaceETA.reset(new ETFESPACE1(M_pressureFESpace->mesh(),
                                           &(M_pressureFESpace->refFE()),
                                             comm));
}

shp<BlockVector>
StokesAssemblerFE::
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

}
