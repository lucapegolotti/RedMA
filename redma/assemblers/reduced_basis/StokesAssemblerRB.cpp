#include "StokesAssemblerRB.hpp"

namespace RedMA
{

StokesAssemblerRB::
StokesAssemblerRB(const DataContainer& data, shp<TreeNode> treeNode) :
  aAssemblerRB(data, treeNode),
  StokesModel(data,treeNode)
{
    this->M_bcManager.reset(new BCManager(data, treeNode));
    M_name = "StokesAssemblerRB";
    M_nComponents = 2;
}

void
StokesAssemblerRB::
setup()
{
    Chrono chrono;
    chrono.start();

    std::string msg = "[";
    msg += this->M_name;
    msg += "] initializing ...";
    printlog(YELLOW, msg, this->M_data.getVerbose());
    initializeFEspaces();
    M_mass = spcast<BlockMatrix>(assembleReducedMass(M_bcManager)); // #1
    M_stiffness = spcast<BlockMatrix>(assembleReducedStiffness(M_bcManager)); // #2
    M_divergence = spcast<BlockMatrix>(assembleReducedDivergence(M_bcManager)); // #3
    assembleFlowRateVectors();
    // assembleFlowRateJacobians();

    setExporter();

    // initializePythonStructures();

    msg = "done, in ";
    msg += std::to_string(chrono.diff());
    msg += " seconds\n";
    printlog(YELLOW, msg, this->M_data.getVerbose());
}

shp<aMatrix>
StokesAssemblerRB::
getMass(const double& time, const shp<aVector>& sol)
{
    return M_reducedMass;
}

shp<aMatrix>
StokesAssemblerRB::
getMassJacobian(const double& time, const shp<aVector>& sol)
{
    shp<BlockMatrix> retMat(new BlockMatrix(2,2));
    return retMat;
}

shp<aVector>
StokesAssemblerRB::
getRightHandSide(const double& time, const shp<aVector>& sol)
{
    shp<BlockMatrix> systemMatrix(new BlockMatrix(this->M_nComponents,
                                                  this->M_nComponents));
    systemMatrix->add(M_reducedStiffness);
    systemMatrix->add(M_reducedDivergence);
    systemMatrix->multiplyByScalar(-1.0);
    shp<aVector> retVec = systemMatrix->multiplyByVector(sol);
    return retVec;
}

shp<aMatrix>
StokesAssemblerRB::
getJacobianRightHandSide(const double& time, const shp<aVector>& sol)
{
    shp<BlockMatrix> retMat(new BlockMatrix(this->M_nComponents, this->M_nComponents));

    retMat->add(M_reducedStiffness);
    retMat->add(M_reducedDivergence);
    retMat->multiplyByScalar(-1.0);

    return retMat;
}

void
StokesAssemblerRB::
initializeFEspaces()
{
    initializeVelocityFESpace(M_comm);
    initializePressureFESpace(M_comm);
}

void
StokesAssemblerRB::
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

std::vector<shp<aMatrix>>
StokesAssemblerRB::
getMatrices() const
{
    std::vector<shp<aMatrix>> retVec;

    retVec.push_back(M_reducedMass);
    retVec.push_back(M_reducedStiffness);
    retVec.push_back(M_reducedDivergence);

    return retVec;
}

shp<aMatrix>
StokesAssemblerRB::
assembleMatrix(const unsigned int& index)
{
    throw new Exception("assembleMatrix not implemented for StokesAssemblerRB");
}

void
StokesAssemblerRB::
applyDirichletBCsMatrix(shp<aMatrix> matrix, double diagCoeff) const
{
    // throw new Exception("Method not implemented for RB");
}

void
StokesAssemblerRB::
applyDirichletBCs(const double& time, shp<aVector> vector) const
{
    // throw new Exception("Method not implemented for RB");
}

void
StokesAssemblerRB::
postProcess(const double& t, const double &dt, const shp<aVector>& sol)
{
    M_bcManager->postProcess();
}

void
StokesAssemblerRB::
apply0DirichletBCs(shp<aVector> vector) const
{
    // throw new Exception("Method not implemented for RB");
}

void
StokesAssemblerRB::
exportSolution(const double& t, const shp<aVector>& sol)
{
    auto solBlck = convert<BlockVector>(sol);
    std::ofstream outfile("rbcoefs/block" + std::to_string(StokesModel::M_treeNode->M_ID) + "t" + std::to_string(t) + ".txt");
    std::string str2write = spcast<DenseVector>(solBlck->block(0))->getString(',') + "\n";
    outfile.write(str2write.c_str(), str2write.size());
    outfile.close();
    unsigned int id = StokesModel::M_treeNode->M_ID;

    *M_velocityExporter = *M_bases->reconstructFEFunction(solBlck->block(0), 0, id);
    *M_pressureExporter = *M_bases->reconstructFEFunction(solBlck->block(1), 1, id);

    bool exportWSS = this->M_data("exporter/export_wss", 1);
    if (exportWSS)
        computeWallShearStress(M_velocityExporter, M_WSSExporter, M_comm);

    // BlockVector<VectorEp> solCopy(2);
    // solCopy.block(0).data() = M_velocityExporter;
    // computeFlowRates(solCopy, true);

    exportNorms(t, M_velocityExporter, M_pressureExporter);

    CoutRedirecter ct;
    ct.redirect();
    M_exporter->postProcess(t);
    printlog(CYAN, ct.restore());
}

shp<aVector>
StokesAssemblerRB::
getZeroVector() const
{
    shp<BlockVector> retVec(new BlockVector(M_nComponents));

    shp<DENSEVECTOR> uComp(new DENSEVECTOR(M_bases->getSizeEnrichedBasis(0)));
    uComp->Scale(0.0);
    shp<DENSEVECTOR> pComp(new DENSEVECTOR(M_bases->getSizeEnrichedBasis(1)));
    pComp->Scale(0.0);

    shp<DenseVector> uCompWrap(new DenseVector());
    uCompWrap->setVector(uComp);

    shp<DenseVector> pCompWrap(new DenseVector());
    pCompWrap->setVector(pComp);

    retVec->setBlock(0,uCompWrap);
    retVec->setBlock(1,pCompWrap);

    return retVec;
}

shp<aVector>
StokesAssemblerRB::
getLifting(const double& time) const
{
    shp<BlockVector> liftingFE(new BlockVector(2));
    // velocity
    shp<DistributedVector> ucomp(new DistributedVector());
    ucomp->setData(shp<VECTOREPETRA>(new VECTOREPETRA(M_velocityFESpace->map(),LifeV::Unique)));
    ucomp->multiplyByScalar(0);

    shp<DistributedVector> pcomp(new DistributedVector());
    pcomp->setData(shp<VECTOREPETRA>(new VECTOREPETRA(M_pressureFESpace->map(),LifeV::Unique)));
    pcomp->multiplyByScalar(0);

    liftingFE->setBlock(0,ucomp);
    liftingFE->setBlock(1,pcomp);

    this->M_bcManager->applyDirichletBCs(time, *liftingFE, this->getFESpaceBCs(),
                                         this->getComponentBCs());

    // shp<aVector> lifting = M_bases->leftProject(liftingFE, StokesModel::M_treeNode->M_ID);
    return liftingFE;
}

void
StokesAssemblerRB::
RBsetup()
{
    if (M_bases == nullptr)
        throw new Exception("RB bases have not been set yet");

    // scale with piola
    unsigned int indexField = 0;
    M_bases->scaleBasisWithPiola(0, StokesModel::M_treeNode->M_ID, [=](shp<VECTOREPETRA> vector)
    {
        shp<BlockVector> vectorWrap(new BlockVector(2));

        shp<DistributedVector> compVec(new DistributedVector());
        compVec->setVector(vector);

        vectorWrap->setBlock(0,compVec);

        applyPiola(vectorWrap, false);
    });

    // // restrict rb matrices
    // if (M_data("rb/online/usemdeim", true))
    // {
    //     std::vector<unsigned int> selectorsU = M_bases->getSelectors(0);
    //     std::vector<unsigned int> selectorsP = M_bases->getSelectors(1);
    //
    //     // this is the case when we do not choose to keep all the vectors in the basis
    //     if (selectorsU.size() > 0)
    //     {
    //         unsigned int Nu = selectorsU.size();
    //         unsigned int Np = selectorsP.size();
    //
    //         // restrict mass
    //         SHP(DENSEMATRIX) restrictedMass(new DENSEMATRIX(Nu,Nu));
    //
    //         unsigned int inew = 0;
    //         for (auto i : selectorsU)
    //         {
    //             unsigned int jnew = 0;
    //             for (auto j : selectorsU)
    //             {
    //                 (*restrictedMass)(inew,jnew) = (*M_mass.block(0,0).data())(i,j);
    //                 jnew++;
    //             }
    //             inew++;
    //         }
    //
    //         M_mass.block(0,0).data() = restrictedMass;
    //
    //         // restrict stiffness
    //         SHP(DENSEMATRIX) restrictedStiffness(new DENSEMATRIX(Nu,Nu));
    //
    //         inew = 0;
    //         for (auto i : selectorsU)
    //         {
    //             unsigned int jnew = 0;
    //             for (auto j : selectorsU)
    //             {
    //                 (*restrictedStiffness)(inew,jnew) = (*M_stiffness.block(0,0).data())(i,j);
    //                 jnew++;
    //             }
    //             inew++;
    //         }
    //
    //         M_stiffness.block(0,0).data() = restrictedStiffness;
    //
    //         // restrict divergence
    //         SHP(DENSEMATRIX) restrictedBT(new DENSEMATRIX(Nu,Np));
    //
    //         inew = 0;
    //         for (auto i : selectorsU)
    //         {
    //             unsigned int jnew = 0;
    //             for (auto j : selectorsP)
    //             {
    //                 (*restrictedBT)(inew,jnew) = (*M_divergence.block(0,1).data())(i,j);
    //                 jnew++;
    //             }
    //             inew++;
    //         }
    //
    //         M_divergence.block(0,1).data() = restrictedBT;
    //
    //         SHP(DENSEMATRIX) restrictedB(new DENSEMATRIX(Np,Nu));
    //
    //         inew = 0;
    //         for (auto i : selectorsP)
    //         {
    //             unsigned int jnew = 0;
    //             for (auto j : selectorsU)
    //             {
    //                 (*restrictedB)(inew,jnew) = (*M_divergence.block(1,0).data())(i,j);
    //                 jnew++;
    //             }
    //             inew++;
    //         }
    //
    //         M_divergence.block(1,0).data() = restrictedB;
    //     }
    // }
    // else
    {
        printlog(YELLOW, "[StokesAssembler] NOT using MDEIM: assembling and projecting matrices\t", M_data.getVerbose());
        Chrono chrono;
        chrono.start();

        // BlockMatrix<MatrixEp> fullMass = assembleReducedMass(nullptr);
        // BlockMatrix<MatrixEp> fullStiffness = assembleReducedStiffness(nullptr);
        // BlockMatrix<MatrixEp> fullDivergence = assembleReducedDivergence(nullptr);

        unsigned int id = StokesModel::M_treeNode->M_ID;

        M_reducedMass.reset(new BlockMatrix(2,2));
        M_reducedStiffness.reset(new BlockMatrix(2,2));
        M_reducedDivergence.reset(new BlockMatrix(2,2));

        M_reducedMass->setBlock(0,0,M_bases->matrixProject(M_mass->block(0,0), 0, 0, id));
        M_reducedStiffness->setBlock(0,0,M_bases->matrixProject(M_stiffness->block(0,0), 0, 0, id));
        M_reducedDivergence->setBlock(0,1,M_bases->matrixProject(M_divergence->block(0,1), 0, 1, id));
        M_reducedDivergence->setBlock(1,0,M_bases->matrixProject(M_divergence->block(1,0), 1, 0, id));

        // M_mass.block(0,0) = M_bases->matrixProject(fullMass.block(0,0), 0, 0, id);
        // M_stiffness.block(0,0) = M_bases->matrixProject(fullStiffness.block(0,0), 0, 0, id);
        // M_divergence.block(0,1) = M_bases->matrixProject(fullDivergence.block(0,1), 0, 1, id);
        // M_divergence.block(1,0) = M_bases->matrixProject(fullDivergence.block(1,0), 1, 0, id);

        std::string msg = "done, in ";
        msg += std::to_string(chrono.diff());
        msg += " seconds\n";
        printlog(YELLOW, msg, this->M_data.getVerbose());
    }
}

shp<RBBases>
StokesAssemblerRB::
getRBBases() const
{
    return M_bases;
}

void
StokesAssemblerRB::
setRBBases(shp<RBBasesManager> rbManager)
{
    // std::string mdeimdir = this->M_data("rb/online/mdeim/directory", "mdeims");
    std::string meshName = StokesModel::M_treeNode->M_block->getMeshName();
    unsigned int dashpos = meshName.find("/");
    unsigned int formatpos = meshName.find(".mesh");
    std::string actualName = meshName.substr(dashpos + 1,
                                             formatpos - dashpos - 1);

    // beware that at this point the rb bases have not been loaded yet
    M_bases = rbManager->getRBBases(actualName);
    M_bases->setFESpace(M_velocityFESpace, 0);
    M_bases->setFESpace(M_pressureFESpace, 1);
}

shp<aVector>
StokesAssemblerRB::
convertFunctionRBtoFEM(shp<aVector> rbSolution) const
{
    auto rbSolutionBlck = convert<BlockVector>(rbSolution);

    shp<BlockVector> retVec(new BlockVector(2));

    unsigned int id = StokesModel::M_treeNode->M_ID;

    if (rbSolutionBlck->block(0)->data())
    {
        shp<DenseVector> comp0 = convert<DenseVector>(rbSolutionBlck->block(0));
        shp<DistributedVector> rec0(new DistributedVector());
        rec0->setVector(M_bases->reconstructFEFunction(comp0, 0, id));
        retVec->setBlock(0,rec0);
    }

    if (rbSolutionBlck->block(1)->data())
    {
        shp<DenseVector> comp1 = convert<DenseVector>(rbSolutionBlck->block(1));
        shp<DistributedVector> rec1(new DistributedVector());
        rec1->setVector(M_bases->reconstructFEFunction(comp1, 1, id));
        retVec->setBlock(1,rec1);
    }

    return retVec;
}

void
StokesAssemblerRB::
applyPiola(shp<aVector> solution, bool inverse)
{
    // std::string msg = "[";
    // msg += this->M_name;
    // msg += "] apply Piola\n";
    // printlog(YELLOW, msg, this->M_data.getVerbose());

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
