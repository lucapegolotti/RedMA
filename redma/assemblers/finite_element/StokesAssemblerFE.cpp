#include "StokesAssemblerFE.hpp"

namespace RedMA
{

StokesAssemblerFE::
StokesAssemblerFE(const DataContainer& data, SHP(TreeNode) treeNode) :
  StokesAssembler(data, treeNode)
{
    this->M_bcManager.reset(new BCManager(data, treeNode));
    M_name = "StokesAssemblerFE";
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

    M_mass = assembleReducedMass(M_bcManager); // #1
    M_stiffness = assembleReducedStiffness(M_bcManager); // #2
    M_divergence = assembleReducedDivergence(M_bcManager); // #3

    assembleFlowRateVectors();
    // assembleFlowRateJacobians();

    setExporter();

    // initializePythonStructures();

    msg = "done, in ";
    msg += std::to_string(chrono.diff());
    msg += " seconds\n";
    printlog(YELLOW, msg, this->M_data.getVerbose());
}

void
StokesAssemblerFE::
postProcess(const double& t, const SHP(BlockVector)& sol)
{
    // shift solutions in multistep method embedded in windkessels
    M_bcManager->postProcess();
}

void
StokesAssemblerFE::
exportSolution(const double& t, const SHP(BlockVector)& sol)
{
    *M_velocityExporter = *static_cast<VECTOREPETRA*>(sol->block(0)->data().get());
    *M_pressureExporter = *static_cast<VECTOREPETRA*>(sol->block(1)->data().get());

    bool exportWSS = this->M_data("exporter/export_wss", 1);
    if (exportWSS)
        computeWallShearStress(M_velocityExporter, M_WSSExporter, M_comm);

    SHP(BlockVector) solCopy(new BlockVector(2));
    solCopy->block(0)->setData(M_velocityExporter);
    StokesAssemblerFE::computeFlowRates(solCopy, true);

    StokesAssembler::exportNorms(t, M_velocityExporter, M_pressureExporter);

    CoutRedirecter ct;
    ct.redirect();
    M_exporter->postProcess(t);
    printlog(CYAN, ct.restore());
}

SHP(BlockMatrix)
StokesAssemblerFE::
getMass(const double& time, const SHP(BlockVector)& sol)
{
    return M_mass;
}

SHP(BlockMatrix)
StokesAssemblerFE::
getMassJacobian(const double& time, const SHP(BlockVector)& sol)
{
    SHP(BlockMatrix) retMat(new BlockMatrix(this->M_nComponents, this->M_nComponents));
    return retMat;
}

SHP(BlockVector)
StokesAssemblerFE::
getRightHandSide(const double& time,
                 const SHP(BlockVector)& sol)
{
    SHP(BlockMatrix) systemMatrix(new BlockMatrix());

    systemMatrix->resize(this->M_nComponents, this->M_nComponents);
    systemMatrix->add(M_stiffness);
    systemMatrix->add(M_divergence);
    systemMatrix->multiplyByScalar(-1.0);

    SHP(BlockVector) retVec;
    retVec.reset(static_cast<BlockVector*>(systemMatrix->multiplyByVector(sol).get()));

    // addNeumannBCs(retVec, time, sol);

    return retVec;
}

SHP(BlockMatrix)
StokesAssemblerFE::
getJacobianRightHandSide(const double& time, const SHP(BlockVector)& sol)
{
    SHP(BlockMatrix) retMat(new BlockMatrix());
    retMat->resize(this->M_nComponents, this->M_nComponents);

    retMat->add(M_stiffness);
    retMat->add(M_divergence);
    retMat->multiplyByScalar(-1.0);

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

SHP(BlockVector)
StokesAssemblerFE::
getZeroVector() const
{
    return buildZeroVector();
}

SHP(BlockVector)
StokesAssemblerFE::
getFELifting(const double& time) const
{
    // BlockVector<FEVECTOR> lifting;
    // lifting.resize(2);
    // lifting.block(0).data().reset(new VECTOREPETRA(M_velocityFESpace->map(),
    //                                                LifeV::Unique));
    // lifting.block(0).data()->zero();
    // lifting.block(1).data().reset(new VECTOREPETRA(M_pressureFESpace->map(),
    //                                                LifeV::Unique));
    // lifting.block(1).data()->zero();
    //
    // this->M_bcManager->applyDirichletBCs(time, lifting, this->getFESpaceBCs(),
    //                                      this->getComponentBCs());
    //
    // return lifting;
}

SHP(BlockVector)
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

void
StokesAssemblerFE::
apply0DirichletBCsMatrix(SHP(BlockMatrix)& matrix, double diagCoeff) const
{
    this->M_bcManager->apply0DirichletMatrix(*matrix, getFESpaceBCs(),
                                             getComponentBCs(), diagCoeff);
}

void
StokesAssemblerFE::
apply0DirichletBCs(SHP(BlockVector)& vector) const
{
    this->M_bcManager->apply0DirichletBCs(*vector, getFESpaceBCs(), getComponentBCs());
}

void
StokesAssemblerFE::
applyDirichletBCs(const double& time, SHP(BlockVector)& vector) const
{
    this->M_bcManager->applyDirichletBCs(time, *vector, getFESpaceBCs(),
                                         getComponentBCs());
}

SHP(FESPACE)
StokesAssemblerFE::
getFEspace(unsigned int index) const
{
    if (index == 0)
        return M_velocityFESpace;
    else if (index == 1)
        return M_pressureFESpace;

    return nullptr;
}

std::vector<SHP(BlockMatrix)>
StokesAssemblerFE::
getMatrices() const
{
    std::vector<SHP(BlockMatrix)> retVec;

    retVec.push_back(M_mass);
    retVec.push_back(M_stiffness);
    retVec.push_back(M_divergence);

    return retVec;
}

SHP(BlockMatrix)
StokesAssemblerFE::
assembleMatrix(const unsigned int& index,
               BlockMDEIMStructure* structure)
{
    if (index == 0)
    {
        return assembleReducedMass(M_bcManager, structure);
    }
    else if (index == 1)
    {
        return assembleReducedStiffness(M_bcManager, structure);
    }
    else if (index == 2)
    {
        return assembleReducedDivergence(M_bcManager, structure);
    }
}

SHP(SparseMatrix)
StokesAssemblerFE::
getNorm(const unsigned int& fieldIndex, bool bcs)
{
    using namespace LifeV;
    using namespace ExpressionAssembly;

    SHP(SparseMatrix) retMat;
    if (fieldIndex == 0)
    {
        // if (!M_massVelocity.data())
        // {
            SHP(MATRIXEPETRA) Nu(new MATRIXEPETRA(M_velocityFESpace->map()));

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
                SHP(BlockMatrix) normWrap(new BlockMatrix());
                normWrap->resize(1,1);
                normWrap->block(0,0)->data() = Nu;

                // note. Applying bcs does not change the norm if Dirichlet bcs are
                // homogeneous (=> lifting) or imposed weakly. Here we impose bcs
                // in order to have the correct conditions in the computation of the
                // supremizers (we have to solve a linear system..)
                apply0DirichletBCsMatrix(normWrap, 1.0);
            }

            M_massVelocity->setData(Nu);
            retMat->setData(Nu);
        // }
        // else
        //     retMat = M_massVelocity;
    }
    else
    {
        // if (!M_massPressure.data())
        // {
            SHP(MATRIXEPETRA) Np(new MATRIXEPETRA(M_pressureFESpace->map()));

            integrate(elements(M_pressureFESpaceETA->mesh()),
                      M_pressureFESpace->qr(),
                      M_pressureFESpaceETA,
                      M_pressureFESpaceETA,
                      phi_i * phi_j
                  ) >> Np;

            Np->globalAssemble();

            M_massPressure->setData(Np);
            retMat->setData(Np);
        // }
        // else
        //     retMat = M_massPressure;
    }

    return retMat;
}

SHP(SparseMatrix)
StokesAssemblerFE::
getConstraintMatrix()
{
    return std::static_pointer_cast<SparseMatrix>(M_divergence->block(0,1));
}

void
StokesAssemblerFE::
setMDEIMs(SHP(MDEIMManager) mdeimManager)
{
   // std::string mdeimdir = this->M_data("rb/online/mdeim/directory", "mdeims");
   //     // std::string meshName = this->M_treeNode->M_block->getMeshName();
   //     // unsigned int dashpos = meshName.find("/");
   //     // unsigned int formatpos = meshName.find(".mesh");
   //     // std::string actualName = meshName.substr(dashpos + 1,
   //     //                                          formatpos - dashpos - 1);
   //     //
   //     // auto mdeims = mdeimManager->getMDEIMS(actualName);
   //     //
   //     // M_mdeimMass = mdeims[0];
   //     // M_mdeimStiffness = mdeims[1];
   //     // M_mdeimDivergence = mdeims[2];
}

void
StokesAssemblerFE::
setExtrapolatedSolution(const SHP(BlockVector)& exSol)
{
    M_extrapolatedSolution->softCopy(exSol);
}

void
StokesAssemblerFE::
applyPiola(SHP(BlockVector) solution, bool inverse)
{
    std::string msg = "[";
    msg += this->M_name;
    msg += "] apply Piola\n";
    printlog(YELLOW, msg, this->M_data.getVerbose());

    using namespace LifeV;
    using namespace ExpressionAssembly;

    auto defAssembler = this->M_defaultAssemblers->
                        getDefaultAssembler(aAssembler::M_treeNode->M_block->getMeshName());

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

    auto velocity = std::static_pointer_cast<VECTOREPETRA>(solution->block(0)->data());

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
            aAssembler::M_treeNode->M_block->matrixInverse(jacobian, transformationMatrix);
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
