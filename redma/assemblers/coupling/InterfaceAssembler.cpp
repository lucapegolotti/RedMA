#include "InterfaceAssembler.hpp"

namespace RedMA
{

shp<LifeV::QuadratureRule>
InterfaceAssembler::
generateQuadratureRule(std::string tag) const
{
    using namespace LifeV;

    shp<QuadratureRule> customRule;
    if (!std::strcmp(tag.c_str(),"STRANG10"))
    {
        // rule taken from https://people.sc.fsu.edu/~jburkardt/datasets/quadrature_rules_tri/quadrature_rules_tri.html
        customRule.reset(new QuadratureRule("STRANG10",TRIANGLE, 3, 7, 0));

        QuadraturePoint p1(0.333333333333333,0.333333333333333,0,
                           -0.149570044467670/2);
        customRule->addPoint(p1);

        QuadraturePoint p2(0.479308067841923,0.260345966079038,0,
                           0.175615257433204/2);
        customRule->addPoint(p2);

        QuadraturePoint p3(0.260345966079038,0.479308067841923,0,
                           0.175615257433204/2);
        customRule->addPoint(p3);

        QuadraturePoint p4(0.260345966079038,0.260345966079038,0,
                           0.175615257433204/2);
        customRule->addPoint(p4);

        QuadraturePoint p5(0.869739794195568,0.065130102902216,0,
                           0.053347235608839/2);
        customRule->addPoint(p5);

        QuadraturePoint p6(0.065130102902216,0.869739794195568,0,
                           0.053347235608839/2);
        customRule->addPoint(p6);

        QuadraturePoint p7(0.065130102902216,0.065130102902216,0,
                           0.053347235608839/2);
        customRule->addPoint(p7);

        QuadraturePoint p8(0.638444188569809,0.312865496004875,0,
                           0.077113760890257/2);
        customRule->addPoint(p8);

        QuadraturePoint p9(0.638444188569809,0.048690315425316,0,
                           0.077113760890257/2);
        customRule->addPoint(p9);

        QuadraturePoint p10(0.312865496004875,0.638444188569809,0,
                            0.077113760890257/2);
        customRule->addPoint(p10);

        QuadraturePoint p11(0.312865496004875,0.048690315425316,0,
                            0.077113760890257/2);
        customRule->addPoint(p11);

        QuadraturePoint p12(0.048690315425316,0.638444188569809,0,
                            0.077113760890257/2);
        customRule->addPoint(p12);

        QuadraturePoint p13(0.048690315425316,0.312865496004875,0,
                            0.077113760890257/2);
        customRule->addPoint(p13);
    }
    else if (!std::strcmp(tag.c_str(),"TOMS612_19"))
    {
        // rule taken from https://people.sc.fsu.edu/~jburkardt/datasets/quadrature_rules_tri/quadrature_rules_tri.html
        customRule.reset(new QuadratureRule("TOMS612_19",TRIANGLE, 3, 9, 0));

        QuadraturePoint p1(0.33333333333333331,0.33333333333333331,0,
                           9.71357962827961025e-2/2);
        customRule->addPoint(p1);

        QuadraturePoint p2(2.06349616025259287e-2,0.48968251919873701,0,
                           3.13347002271398278e-2/2);
        customRule->addPoint(p2);

        QuadraturePoint p3(0.48968251919873701,2.06349616025259287e-2,0,
                           3.13347002271398278e-2/2);
        customRule->addPoint(p3);

        QuadraturePoint p4(0.48968251919873701,0.48968251919873701,0,
                           3.13347002271398278e-2/2);
        customRule->addPoint(p4);

        QuadraturePoint p5(0.12582081701412900,0.43708959149293553,0,
                           7.78275410047754301e-2/2);
        customRule->addPoint(p5);

        QuadraturePoint p6(0.43708959149293553,0.12582081701412900,0,
                           7.78275410047754301e-2/2);
        customRule->addPoint(p6);

        QuadraturePoint p7(0.43708959149293553,0.43708959149293553,0,
                           7.78275410047754301e-2/2);
        customRule->addPoint(p7);

        QuadraturePoint p8(0.62359292876193562,0.18820353561903219,0,
                           7.96477389272090969e-2/2);
        customRule->addPoint(p8);

        QuadraturePoint p9(0.18820353561903219,0.62359292876193562,0,
                           7.96477389272090969e-2/2);
        customRule->addPoint(p9);

        QuadraturePoint p10(0.18820353561903219,0.18820353561903219,0,
                            7.96477389272090969e-2/2);
        customRule->addPoint(p10);

        QuadraturePoint p11(0.91054097321109406,4.47295133944529688e-2,0,
                            2.55776756586981006e-2/2);
        customRule->addPoint(p11);

        QuadraturePoint p12(4.47295133944529688e-2,0.91054097321109406,0,
                            2.55776756586981006e-2/2);
        customRule->addPoint(p12);

        QuadraturePoint p13(4.47295133944529688e-2,4.47295133944529688e-2,0,
                            2.55776756586981006e-2/2);
        customRule->addPoint(p13);

        QuadraturePoint p14(0.74119859878449801,3.68384120547362581e-2,0,
                            4.32835393772893970e-2/2);
        customRule->addPoint(p14);

        QuadraturePoint p15(0.74119859878449801,0.22196298916076573,0,
                            4.32835393772893970e-2/2);
        customRule->addPoint(p15);

        QuadraturePoint p16(3.6838412054736258e-2,0.74119859878449801,0,
                            4.32835393772893970e-2/2);
        customRule->addPoint(p1);

        QuadraturePoint p17(3.68384120547362581e-2,0.22196298916076573,0,
                            4.32835393772893970e-2/2);
        customRule->addPoint(p17);

        QuadraturePoint p18(0.22196298916076573,0.74119859878449801,0,
                            4.32835393772893970e-2/2);
        customRule->addPoint(p18);

        QuadraturePoint p19(0.22196298916076573,3.68384120547362581e-2,0,
                            4.32835393772893970e-2/2);
        customRule->addPoint(p19);
    }
    else
    {
        throw new Exception("Quadrature rule " + tag + " not implemented!");
    }
    return customRule;
}

Interface::
Interface() :
  M_inletIndex(0)
{
}

Interface::
Interface(shp<AssemblerType> assemblerFather,
          const int& indexFather,
          shp<AssemblerType> assemblerChild,
          const int& indexChild,
          const unsigned int& interfaceID) :
  M_assemblerFather(assemblerFather),
  M_indexFather(indexFather),
  M_assemblerChild(assemblerChild),
  M_indexChild(indexChild),
  M_ID(interfaceID),
  M_inletIndex(0)
{
}

InterfaceAssembler::
InterfaceAssembler(const DataContainer& data,
                   const bool& addNoSlipBC) :
  M_data(data),
  M_isInlet(false),
  M_addNoSlipBC(addNoSlipBC)
{
}

InterfaceAssembler::
InterfaceAssembler(const DataContainer& data,
                   const Interface& interface,
                   const bool& addNoSlipBC) :
  M_data(data),
  M_interface(interface),
  M_isInlet(false),
  M_addNoSlipBC(addNoSlipBC)
{
    if (M_interface.M_indexFather == -1 || M_interface.M_indexChild == -1)
        M_isInlet = true;
    setup();
}

void
InterfaceAssembler::
setup()
{
    Chrono chrono;
    chrono.start();

    printlog(YELLOW, "[InterfaceAssembler] initialize interface"
                     " assembler ...", M_data.getVerbose());

    M_stabilizationCoupling = M_data("coupling/stab_coefficient", 0.0);
    // we set size to zero because it will be reset in the following method
    M_fatherBT.reset(new BlockMatrix(0,0));
    M_fatherB.reset(new BlockMatrix(0,0));
    M_childBT.reset(new BlockMatrix(0,0));
    M_childB.reset(new BlockMatrix(0,0));
    buildCouplingMatrices();

    std::string msg = "done, in ";
    msg += std::to_string(chrono.diff());
    msg += " seconds\n";
    printlog(YELLOW, msg, M_data.getVerbose());
}

std::vector<shp<DistributedVector>>
InterfaceAssembler::
buildCouplingVectors(shp<BasisFunctionFunctor> bfs,
                     const GeometricFace& face,
                     shp<aAssembler> assembler) const
{
    using namespace LifeV;
    using namespace ExpressionAssembly;

    QuadratureBoundary boundaryQuadRule(buildTetraBDQR
                                        (*generateQuadratureRule("STRANG10")));

    unsigned int nBasisFunctions = bfs->getNumBasisFunctions();

    std::vector<shp<DistributedVector>> couplingVectors(3 * nBasisFunctions);

    shp<ETFESPACE3> etfespace = assembler->getETFESpaceCoupling();
    MAPEPETRA map = etfespace->map();
    shp<MESH> mesh = etfespace->mesh();
    unsigned int faceFlag = face.M_flag;

    unsigned int count = 0;
    LifeV::VectorSmall<3> versor;

    // this must be changed for scalar problems (e.g. laplacian)
    for (unsigned int dim = 0; dim < 3; dim++)
    {
        versor *= 0.;
        versor[dim] = 1.;

        for (unsigned int i = 0; i < nBasisFunctions; i++)
        {
            shp<VECTOREPETRA> currentMode(new VECTOREPETRA(map, LifeV::Repeated));

            bfs->setIndex(i);
            couplingVectors[count].reset(new DistributedVector());

            integrate(boundary(mesh, faceFlag),
                      boundaryQuadRule,
                      etfespace,
                      eval(bfs, X) * dot(versor, phi_i)
                  ) >> currentMode;

            couplingVectors[count]->setData(currentMode);

            count++;
        }
    }

    return couplingVectors;
}

void
InterfaceAssembler::
buildCouplingMatrices(shp<AssemblerType> assembler,
                      const GeometricFace& face,
                      shp<BlockMatrix> matrixT,
                      shp<BlockMatrix> matrix,
                      const bool isFather)
{
    shp<BasisFunctionFunctor> bfs;

    matrixT->resize(assembler->getNumComponents(),1);
    matrix->resize(1, assembler->getNumComponents());

    // 1) add disk basis functions
    bfs = BasisFunctionFactory(M_data.getDatafile(), face, M_isInlet, false);

    std::vector<shp<DistributedVector>> couplingVectors;
    couplingVectors = buildCouplingVectors(bfs, face, assembler);

    // 2) add extra basis function on the disk ring, if necessary
    const double h = (double) LifeV::MeshUtility::MeshStatistics::computeSize(*(assembler->getFEspace(0)->mesh())).meanH;

    bfs = BasisFunctionFactory(M_data.getDatafile(), face, M_isInlet, true, h);

    std::vector<shp<DistributedVector>> ringCouplingVectors;
    ringCouplingVectors = buildCouplingVectors(bfs, face, assembler);

    couplingVectors.insert(couplingVectors.end(),
                           ringCouplingVectors.begin(), ringCouplingVectors.end());

    bool strongRingCoupling = M_data("coupling/strong_ring", false);
    if (!M_isInlet && strongRingCoupling)
    {
        M_mapLagrange = spcast<VECTOREPETRA>(couplingVectors[0]->data())->mapPtr();

        std::vector<shp<DistributedVector>> strongRingCouplingVectors;
        strongRingCouplingVectors = (isFather) ? buildStrongRingCouplingVectors().first :
                                    buildStrongRingCouplingVectors().second;

        couplingVectors.insert(couplingVectors.end(),
                               strongRingCouplingVectors.begin(), strongRingCouplingVectors.end());
    }

    // assembling coupling matrices
    shp<SparseMatrix> couplingMatrix(new SparseMatrix(couplingVectors));
    matrixT->setBlock(0,0,couplingMatrix);
    matrix->setBlock(0,0,couplingMatrix->transpose());

    // imposing BCs on the transpose coupling matrix only
    assembler->getBCManager()->apply0DirichletMatrix(*matrixT,
                                                     assembler->getFESpaceBCs(),
                                                     assembler->getComponentBCs(),
                                                     0.0, !(M_addNoSlipBC));
}

void
InterfaceAssembler::
addContributionRhs(const double& time,
                   shp<BlockVector> rhs,
                   shp<BlockVector> sol,
                   const unsigned int& nPrimalBlocks)
{
    unsigned int fatherID = M_interface.M_indexFather;
    unsigned int childID = M_interface.M_indexChild;
    unsigned int interfaceID = M_interface.M_ID;
    shp<aAssembler> assemblerFather = M_interface.M_assemblerFather;
    shp<aAssembler> assemblerChild = M_interface.M_assemblerChild;

    shp<aVector> tempResFather;
    // we have (-1) because we are solving H un+1 = F(.) and coupling is in F
    tempResFather = M_fatherBT->multiplyByVector(sol->block(nPrimalBlocks + interfaceID));
    tempResFather->multiplyByScalar(-1.0);
    rhs->block(fatherID)->add(tempResFather);

    auto tempResChild = M_childBT->multiplyByVector(sol->block(nPrimalBlocks + interfaceID));
    tempResChild->multiplyByScalar(-1.0);
    rhs->block(childID)->add(tempResChild);

    tempResFather = M_fatherB->multiplyByVector(sol->block(fatherID));
    tempResFather->multiplyByScalar(-1.0);
    if (rhs->block(nPrimalBlocks + interfaceID)->isZero())
        rhs->setBlock(nPrimalBlocks + interfaceID,tempResFather);
    else
        rhs->block(nPrimalBlocks + interfaceID)->add(tempResFather);

    tempResChild = M_childB->multiplyByVector(sol->block(childID));
    tempResChild->multiplyByScalar(-1.0);
    if (rhs->block(nPrimalBlocks + interfaceID)->isZero())
        rhs->setBlock(nPrimalBlocks + interfaceID,tempResChild);
    else
        rhs->block(nPrimalBlocks + interfaceID)->add(tempResChild);
}

double
InterfaceAssembler::
checkStabilizationTerm(const shp<BlockVector>& sol,
                       const unsigned int& nPrimalBlocks)
{
    // if (M_stabilizationCoupling > THRESHOLDSTAB &&
    //     M_interface.M_assemblerFather && M_interface.M_assemblerChild)
    // {
    //     unsigned int fatherID = M_interface.M_indexFather;
    //     unsigned int childID = M_interface.M_indexChild;
    //     unsigned int interfaceID = M_interface.M_ID;
    //
    //     BlockVector<BlockVector<InVectorType>> res;
    //     res.resize(1);
    //
    //     res.block(0) -= (M_stabFather * sol.block(fatherID)) * 0.5;
    //     res.block(0) -= (M_stabChild * sol.block(childID)) * 0.5;
    //
    //     // std::cout << "--------" << std::endl << std::flush;
    //     // res.block(0).block(0).data()->showMe();
    //     // std::cout << "++++++++" << std::endl << std::flush;
    //     // sol.block(nPrimalBlocks + interfaceID).block(0).data()->showMe();
    //     //
    //     // std::cout << "stab term stress " << res.norm2() << std::endl << std::flush;
    //     // std::cout << "stab term lagrange " << sol.block(nPrimalBlocks + interfaceID).norm2() << std::endl << std::flush;
    //
    //     res.block(0) -= sol.block(nPrimalBlocks + interfaceID);
    //
    //     std::string msg = "[InterfaceAssembler] interface ID = ";
    //     msg += std::to_string(interfaceID);
    //     msg += ", stab term norm = ";
    //     msg += std::to_string(res.norm2());
    //     msg += "\n";
    //     printlog(MAGENTA, msg, M_data.getVerbose());
    // }
}

void
InterfaceAssembler::
addContributionJacobianRhs(const double& time,
                           shp<BlockMatrix> jac,
                           shp<BlockVector> sol,
                           const unsigned int& nPrimalBlocks)
{
    unsigned int fatherID = M_interface.M_indexFather;
    unsigned int childID = M_interface.M_indexChild;
    unsigned int interfaceID = M_interface.M_ID;

    // hard copy, otherwise we flip the sign of the matrices every time this function is called
    jac->setBlock(fatherID, nPrimalBlocks + interfaceID, shp<BlockMatrix>(M_fatherBT->clone()));
    jac->setBlock(childID, nPrimalBlocks + interfaceID, shp<BlockMatrix>(M_childBT->clone()));
    jac->setBlock(nPrimalBlocks + interfaceID, fatherID, shp<BlockMatrix>(M_fatherB->clone()));
    jac->setBlock(nPrimalBlocks + interfaceID, childID, shp<BlockMatrix>(M_childB->clone()));

    jac->block(fatherID, nPrimalBlocks + interfaceID)->multiplyByScalar(-1);
    jac->block(childID,  nPrimalBlocks + interfaceID)->multiplyByScalar(-1);
    jac->block(nPrimalBlocks + interfaceID, fatherID)->multiplyByScalar(-1);
    jac->block(nPrimalBlocks + interfaceID,  childID)->multiplyByScalar(-1);
}

shp<BlockVector>
InterfaceAssembler::
getZeroVector() const
{
    shp<BlockVector> retVector(new BlockVector(1));

    shp<VECTOREPETRA> zeroVec(new VECTOREPETRA(*M_mapLagrange, LifeV::Unique));
    zeroVec->zero();

    retVector->setBlock(0,wrap(zeroVec));

    return retVector;
}

void
InterfaceAssembler::
buildCouplingMatrices()
{
    unsigned int indexOutlet = M_interface.M_indexOutlet;

    auto asFather = M_interface.M_assemblerFather;
    if (asFather)
    {
        GeometricFace outlet = asFather->getTreeNode()->M_block->getOutlet(indexOutlet);

        buildCouplingMatrices(asFather, outlet, M_fatherBT, M_fatherB, true);

        /*if (M_stabilizationCoupling > THRESHOLDSTAB)
            buildStabilizationMatrix(asFather, outlet, M_stabFather);*/

        M_mapLagrange = spcast<MATRIXEPETRA>(M_fatherBT->block(0,0)->data())->domainMapPtr();

        if (!std::strcmp(asFather->getTreeNode()->M_block->getDiscretizationMethod().c_str(),"rb"))
        {
            M_fatherBT = asFather->getRBBases()->leftProject(M_fatherBT, asFather->ID());
            M_fatherB = asFather->getRBBases()->rightProject(M_fatherB, asFather->ID());
        }
    }

    auto asChild = M_interface.M_assemblerChild;
    if (asChild)
    {
        GeometricFace inlet = asChild->getTreeNode()->M_block->getInlet(M_interface.M_inletIndex);
        // I invert the normal of the face such that it is the same as the outlet
        inlet.M_normal *= (-1.);

        buildCouplingMatrices(asChild, inlet, M_childBT, M_childB, false);

        /*if (M_stabilizationCoupling > THRESHOLDSTAB)
            buildStabilizationMatrix(asChild, inlet, M_stabChild);*/

        M_childBT->multiplyByScalar(-1.);
        M_childB->multiplyByScalar(-1.);

        M_childBTfe = M_childBT;
        M_childBfe = M_childB;

        M_mapLagrange = spcast<MATRIXEPETRA>(M_childBT->block(0,0)->data())->domainMapPtr();

        if (!std::strcmp(asChild->getTreeNode()->M_block->getDiscretizationMethod().c_str(),"rb"))
        {
            M_childBT = asChild->getRBBases()->leftProject(M_childBT, asChild->ID());
            M_childB = asChild->getRBBases()->rightProject(M_childB, asChild->ID());
        }
    }
}

void
InterfaceAssembler::
buildStabilizationMatrix(shp<AssemblerType> assembler,
                         const GeometricFace& face,
                         shp<BlockMatrix> matrix)
{
    // BlockMatrix<MatrixEp> stab;
    // stab.resize(1, assembler->getNumComponents());
    //
    // shp<BasisFunctionFunctor> bfs;
    // bfs =  BasisFunctionFactory(M_data.getDatafile(), face);
    //
    // std::vector<VectorEp> stabVectorsVelocity;
    // stabVectorsVelocity = buildStabilizationVectorsVelocity(bfs, face, assembler);
    //
    // std::vector<VectorEp> stabVectorsPressure;
    // stabVectorsPressure = buildStabilizationVectorsPressure(bfs, face, assembler);
    //
    // stab.block(0,0).shallowCopy(MatrixEp(stabVectorsVelocity).transpose());
    // stab.block(0,1).shallowCopy(MatrixEp(stabVectorsPressure).transpose());
    // matrix.shallowCopy(stab);
    //
    // std::vector<VectorEp> stabVectorsLagrange = buildStabilizationVectorsLagrange();
    // M_identity.resize(1,1);
    //
    // M_identity.block(0,0).shallowCopy(MatrixEp(stabVectorsLagrange));
}

/*InterfaceAssembler::InterfacePtrType
InterfaceAssembler::
buildRingInterfaceMap()
{
    auto asFather = M_interface.M_assemblerFather;
    auto asChild = M_interface.M_assemblerChild;

    shp<FESPACE> FESpaceVFather = asFather->getFEspace(0);
    shp<FESPACE> FESpaceVChild = asChild->getFEspace(0);

    unsigned int indexOutlet = M_interface.M_indexOutlet;
    unsigned int fatherFlag = asFather->getTreeNode()->M_block->getOutlet(indexOutlet).M_flag;
    unsigned int childFlag = asChild->getTreeNode()->M_block->getInlet().M_flag;

    const double tol = 1e-1;

    *//*const double tol = std::max((double) LifeV::MeshUtility::MeshStatistics::computeSize(*(FESpaceVFather->mesh())).maxH,
                                (double) LifeV::MeshUtility::MeshStatistics::computeSize(*(FESpaceVChild->mesh())).maxH);*//*

    InterfacePtrType ifPtr(new InterfaceType());

    ifPtr->setup (FESpaceVFather->refFE(),
                  FESpaceVFather->dof(),
                  FESpaceVChild->refFE(),
                  FESpaceVChild->dof() );
    ifPtr->update (*(FESpaceVFather->mesh()), fatherFlag,
                   *(FESpaceVChild->mesh()), childFlag,
                   tol);

    *//*if (M_data.getVerbose())
        ifPtr->showMe();*//*

    return ifPtr;
}*/

/*std::vector<LifeV::ID>
InterfaceAssembler::
identifyValidRingDOFs()
{
    auto asFather = M_interface.M_assemblerFather;
    if (!asFather)
        throw new Exception("Father assembler is not defined!");

    shp<FESPACE> FESpaceVFather = asFather->getFEspace(0);
    unsigned int indexOutlet = M_interface.M_indexOutlet;
    unsigned int fatherFlag = asFather->getTreeNode()->M_block->getOutlet(indexOutlet).M_ringFlag;
    shp<VECTOREPETRA> ringsIndicator = asFather->getBCManager()->computeRingsIndicator(FESpaceVFather,
                                                                                       fatherFlag, false);

    // Building the father->child DOF map at the interface ring
    InterfacePtrType ifPtr = this->buildRingInterfaceMap();

    std::list<std::pair<LifeV::ID, LifeV::ID>> facetToFacetConnections = ifPtr->connectedFacetMap();
    std::vector<LifeV::ID> localToGlobalMapOnFatherFacets;
    for (auto i = facetToFacetConnections.begin(); i != facetToFacetConnections.end(); ++i)
    {
        std::vector<LifeV::ID> tmpLocalToGlobalMap = FESpaceVFather->dofPtr()->localToGlobalMapOnBdFacet(i->first);

        localToGlobalMapOnFatherFacets.insert(localToGlobalMapOnFatherFacets.end(),
                                              tmpLocalToGlobalMap.begin(), tmpLocalToGlobalMap.end());
    }

    std::sort(localToGlobalMapOnFatherFacets.begin(), localToGlobalMapOnFatherFacets.end());
    localToGlobalMapOnFatherFacets.erase(std::unique(localToGlobalMapOnFatherFacets.begin(),
                                                     localToGlobalMapOnFatherFacets.end()),
                                         localToGlobalMapOnFatherFacets.end());

    std::vector<LifeV::Int> ringIndicesFather;
    unsigned int count = 0;
    for (auto const& fatherID : localToGlobalMapOnFatherFacets)
    {
        if ((ifPtr->isMyInterfaceDof(fatherID) &&
             (std::abs((*ringsIndicator)[fatherID] - fatherFlag) <= 1e-10)))
            ringIndicesFather.push_back(count);
        count++;
    }

    std::vector<LifeV::ID> tmpLocalToGlobalMapOnFatherFacets;
    for (const LifeV::Int ringIndexFather : ringIndicesFather)
        tmpLocalToGlobalMapOnFatherFacets.push_back(localToGlobalMapOnFatherFacets[ringIndexFather]);

    int dofsPerRing = M_data("coupling/strong_dofs_per_ring", 4);
    if (dofsPerRing != -1)
    {
        if ((dofsPerRing & 3) != 0)
            throw new Exception("The number of DOFs to strongly couple at each ring must be a multiple of 4! "
                                "Alternatively, if it equals -1, all ring DOFs are considered when performing the "
                                "strong coupling!");
        else if (dofsPerRing > tmpLocalToGlobalMapOnFatherFacets.size())
            throw new Exception("The number of ring DOFs at which performing the strong coupling exceeds the "
                                "total number of ring DOFs!");
    }

    unsigned int num = (dofsPerRing-4) / 4;
    unsigned int dim = (tmpLocalToGlobalMapOnFatherFacets.size()-4) / 4;

    std::vector<LifeV::ID> validLocalToGlobalMapOnFatherFacets;
    if (dofsPerRing == -1)
        validLocalToGlobalMapOnFatherFacets = tmpLocalToGlobalMapOnFatherFacets;
    else
    {
        validLocalToGlobalMapOnFatherFacets.insert(validLocalToGlobalMapOnFatherFacets.end(),
                                                   tmpLocalToGlobalMapOnFatherFacets.begin(),
                                                   tmpLocalToGlobalMapOnFatherFacets.begin() + 4);

        unsigned int index;
        for (unsigned int quadr = 0; quadr < 4; ++quadr)
            for (unsigned int cnt = 1; cnt <= num; ++cnt)
            {
                index = 3 + (dim * quadr) + (cnt * dim / (num+1)) ;
                validLocalToGlobalMapOnFatherFacets.push_back(tmpLocalToGlobalMapOnFatherFacets[index]);
            }
    }

    return validLocalToGlobalMapOnFatherFacets;
}*/

double
InterfaceAssembler::evaluate_RBF(double x, double y, double z, double xc, double yc, double zc, double R)
{
    double D = std::sqrt(std::pow(x-xc, 2) + std::pow(y-yc, 2) + std::pow(z-zc, 2));
    double val = std::pow(1 - D/R, 4) * (4 * D/R + 1);
    return val;
}

void
InterfaceAssembler::findRingPointsCoordinates()
{
    auto asFather = M_interface.M_assemblerFather;
    auto asChild = M_interface.M_assemblerChild;

    if (asFather && asChild)
    {
        shp<FESPACE > FESpaceVFather = asFather->getFEspace(0);
        shp<FESPACE > FESpaceVChild = asChild->getFEspace(0);

        int numTotalDofFather = FESpaceVFather->dof().numTotalDof();
        unsigned int fatherFlag = asFather->getTreeNode()->M_block->getOutlet(M_interface.M_indexOutlet).M_ringFlag;
        for (unsigned int i = 0; i < numTotalDofFather; i++) {
            if (FESpaceVFather->mesh()->point(i).markerID() == fatherFlag) {
                std::vector<double> coords(3);
                coords[0] = FESpaceVFather->mesh()->point(i).x();
                coords[1] = FESpaceVFather->mesh()->point(i).y();
                coords[2] = FESpaceVFather->mesh()->point(i).z();

                M_fatherRingPoints[FESpaceVFather->mesh()->point(i).id()] = coords;
            }
        }

        int numTotalDofChild = FESpaceVFather->dof().numTotalDof();
        unsigned int childFlag = asChild->getTreeNode()->M_block->getInlet().M_ringFlag;
        for (unsigned int i = 0; i < numTotalDofChild; i++) {
            if (FESpaceVChild->mesh()->point(i).markerID() == childFlag) {
                std::vector<double> coords(3);
                coords[0] = FESpaceVChild->mesh()->point(i).x();
                coords[1] = FESpaceVChild->mesh()->point(i).y();
                coords[2] = FESpaceVChild->mesh()->point(i).z();

                M_childRingPoints[FESpaceVChild->mesh()->point(i).id()] = coords;
            }
        }
    }
}

void
InterfaceAssembler::
buildRingDOFsMap()
{
    for (auto itF = M_fatherRingPoints.begin(); itF != M_fatherRingPoints.end(); ++itF)
    {
        double m = std::numeric_limits<double>::max();
        double d;
        LifeV::ID id = 0;

        for (auto itC= M_childRingPoints.begin(); itC != M_childRingPoints.end(); ++itC)
        {
            d = std::sqrt(std::pow((itF->second[0] - itC->second[0]), 2) +
                          std::pow((itF->second[1] - itC->second[1]), 2) +
                          std::pow((itF->second[2] - itC->second[2]), 2));
            if (d < m) {
                m = d;
                id = itC->first;
            }
        }

        M_ringDOFsMap[itF->first] = id;
    }
}


std::pair<std::vector<shp<DistributedVector>>, std::vector<shp<DistributedVector>>>
InterfaceAssembler::
buildStrongRingCouplingVectors()
{
    auto asFather = M_interface.M_assemblerFather;
    auto asChild = M_interface.M_assemblerChild;

    std::vector<shp<DistributedVector>> SCVectorsFather;
    std::vector<shp<DistributedVector>> SCVectorsChild;

    if (asFather && asChild)
    {
        shp<FESPACE> FESpaceVFather = asFather->getFEspace(0);
        shp<FESPACE> FESpaceVChild = asChild->getFEspace(0);

        // Building the father->child DOF map at the interface ring
        // InterfacePtrType ifPtr = this->buildRingInterfaceMap();

        unsigned int fathernDOF = FESpaceVFather->dofPtr()->numTotalDof();
        unsigned int childnDOF = FESpaceVChild->dofPtr()->numTotalDof();

        const double R = 1.5 * std::max((double) LifeV::MeshUtility::MeshStatistics::computeSize(*(FESpaceVFather->mesh())).maxH,
                                        (double) LifeV::MeshUtility::MeshStatistics::computeSize(*(FESpaceVChild->mesh())).maxH);

        // std::vector<LifeV::ID> validRingDOFs = this->identifyValidRingDOFs();

        this->findRingPointsCoordinates();
        this->buildRingDOFsMap();

        std::vector<LifeV::UInt> indices(2);
        shp<BlockVector> vecFather;
        shp<BlockVector> vecChild;

        LifeV::ID fatherID;
        LifeV::ID childID;
        //for (auto const& fatherID : validRingDOFs)
        for (auto itF = M_fatherRingPoints.begin(); itF != M_fatherRingPoints.end(); ++itF)
        {
            fatherID = itF->first;
            childID = M_ringDOFsMap[fatherID];

            // childID = ifPtr->getInterfaceDof(fatherID);

            /*std::cout << fatherID << " " << childID << std::endl << std::flush;
            std::cout << M_fatherRingPoints[fatherID][0] << " " << M_fatherRingPoints[fatherID][1] << " " << M_fatherRingPoints[fatherID][2] << std::endl;
            std::cout << M_childRingPoints[childID][0] << " " << M_childRingPoints[childID][1] << " " << M_childRingPoints[childID][2] << std::endl;*/

            double val = evaluate_RBF(M_childRingPoints[childID][0],
                                      M_childRingPoints[childID][1],
                                      M_childRingPoints[childID][2],
                                      M_fatherRingPoints[fatherID][0],
                                      M_fatherRingPoints[fatherID][1],
                                      M_fatherRingPoints[fatherID][2],
                                      R);

            for (unsigned int nDim = 0; nDim < 3; ++nDim)
            {
                vecFather = this->getZeroVector();
                vecChild = this->getZeroVector();

                indices[0] = fatherID + nDim * fathernDOF;
                indices[1] = childID + nDim * childnDOF;

                spcast<VECTOREPETRA>(vecFather->block(0)->data())->setCoefficient(indices[0],
                                                                                  val);
                spcast<VECTOREPETRA>(vecChild->block(0)->data())->setCoefficient(indices[1],
                                                                                 1.0);

                SCVectorsFather.push_back(spcast<DistributedVector>(vecFather->block(0)));
                SCVectorsChild.push_back(spcast<DistributedVector>(vecChild->block(0)));
            }
        }
    }

    return std::make_pair(SCVectorsFather, SCVectorsChild);
}

void
InterfaceAssembler::
buildMapLagrange(shp<BasisFunctionFunctor> bfs)
{
    // NOT IMPLEMENTED
}

}  // namespace RedMA
