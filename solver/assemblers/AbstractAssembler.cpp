#include <AbstractAssembler.hpp>

namespace RedMA
{

AbstractAssembler::
AbstractAssembler(const GetPot& datafile, commPtr_Type comm,
                  const TreeNodePtr& treeNode, bool verbose) :
  M_datafile(datafile),
  M_comm(comm),
  M_coupler(comm, verbose),
  M_treeNode(treeNode),
  M_verbose(verbose)
{
}

void
AbstractAssembler::
addPrimalMaps(MapEpetraPtr& globalMap,
              std::vector<MapEpetraPtr>& maps,
              std::vector<unsigned int>& dimensions)
{
    for (MapVectorSTD::iterator it = M_primalMaps.begin();
         it != M_primalMaps.end(); it++)
    {
        *globalMap += *(*it);
        maps.push_back(*it);
        dimensions.push_back((*it)->map(LifeV::Unique)->NumGlobalElements());
    }
}

std::vector<AbstractAssembler::MapEpetraPtr>
AbstractAssembler::
getPrimalMapVector()
{
    return M_primalMaps;
}

std::vector<AbstractAssembler::MapEpetraPtr>
AbstractAssembler::
getDualMapVector()
{
    return M_dualMaps;
}

void
AbstractAssembler::
setTimeIntegrationOrder(unsigned int order)
{
    M_timeIntegrationOrder = order;
}


AbstractAssembler::MatrixPtr
AbstractAssembler::
assembleBoundaryMatrix(GeometricFace face)
{
    using namespace LifeV;
    using namespace ExpressionAssembly;

    // rule taken from https://people.sc.fsu.edu/~jburkardt/datasets/quadrature_rules_tri/quadrature_rules_tri.html
    QuadratureRule customRule("STRANG10",TRIANGLE, 3, 7, 0);

    QuadraturePoint p1(0.333333333333333,0.333333333333333,0,
                       -0.149570044467670/2);
    customRule.addPoint(p1);

    QuadraturePoint p2(0.479308067841923,0.260345966079038,0,
                       0.175615257433204/2);
    customRule.addPoint(p2);

    QuadraturePoint p3(0.260345966079038,0.479308067841923,0,
                       0.175615257433204/2);
    customRule.addPoint(p3);

    QuadraturePoint p4(0.260345966079038,0.260345966079038,0,
                       0.175615257433204/2);
    customRule.addPoint(p4);

    QuadraturePoint p5(0.869739794195568,0.065130102902216,0,
                       0.053347235608839/2);
    customRule.addPoint(p5);

    QuadraturePoint p6(0.065130102902216,0.869739794195568,0,
                       0.053347235608839/2);
    customRule.addPoint(p6);

    QuadraturePoint p7(0.065130102902216,0.065130102902216,0,
                       0.053347235608839/2);
    customRule.addPoint(p7);

    QuadraturePoint p8(0.638444188569809,0.312865496004875,0,
                       0.077113760890257/2);
    customRule.addPoint(p8);

    QuadraturePoint p9(0.638444188569809,0.048690315425316,0,
                       0.077113760890257/2);
    customRule.addPoint(p9);

    QuadraturePoint p10(0.312865496004875,0.638444188569809,0,
                       0.077113760890257/2);
    customRule.addPoint(p10);

    QuadraturePoint p11(0.312865496004875,0.048690315425316,0,
                       0.077113760890257/2);
    customRule.addPoint(p11);

    QuadraturePoint p12(0.048690315425316,0.638444188569809,0,
                       0.077113760890257/2);
    customRule.addPoint(p12);

    QuadraturePoint p13(0.048690315425316,0.312865496004875,0,
                       0.077113760890257/2);
    customRule.addPoint(p13);

    // QuadratureBoundary boundaryQuadRule(buildTetraBDQR(quadRuleTria7pt));
    QuadratureBoundary boundaryQuadRule(buildTetraBDQR(customRule));

    unsigned int faceFlag = face.M_flag;
    MeshPtr mesh = M_couplingFESpaceETA->mesh();

    MapEpetra couplingMap = M_couplingFESpace->map();
    MatrixPtr boundaryMassMatrix(new Matrix(couplingMap));

    // assemble boundary mass matrix to orthonormalize w.r.t. L2 product
    integrate(boundary(mesh, faceFlag),
              boundaryQuadRule,
              M_couplingFESpaceETA,
              M_couplingFESpaceETA,
              dot(phi_i, phi_j)
              ) >> boundaryMassMatrix;
    boundaryMassMatrix->globalAssemble();

    return boundaryMassMatrix;
}

std::shared_ptr<BasisFunctionFunctor>
AbstractAssembler::
castBasisFunction(GeometricFace inlet)
{
    return BasisFunctionFactory(M_datafile, inlet);
}

void
AbstractAssembler::
assembleCouplingMatricesInterpolation(AbstractAssembler& child,
                                      const unsigned int& indexOutlet,
                                      const unsigned int& interfaceIndex,
                                      std::shared_ptr<BasisFunctionFunctor> basisFunction,
                                      MapEpetraPtr& globalMap,
                                      std::vector<MapEpetraPtr>& maps,
                                      std::vector<unsigned int>& dimensions)
{
    AbstractAssembler* mainDomain;
    AbstractAssembler* otherDomain;
    // choose side with smaller h for interpolation
    double hFatherLocal =
        LifeV::MeshUtility::MeshStatistics::computeSize(
                                        *M_couplingFESpace->mesh()).maxH;
    double hChildLocal =
        LifeV::MeshUtility::MeshStatistics::computeSize(
                                  *child.M_couplingFESpace->mesh()).maxH;

    // find global hmax
    double hFather;
    double hChild;

    M_comm->MaxAll(&hFatherLocal, &hFather, 1);
    M_comm->MaxAll(&hChildLocal, &hChild, 1);

    GeometricFace inlet = child.M_treeNode->M_block->getInlet();
    GeometricFace outlet = M_treeNode->M_block->getOutlet(indexOutlet);

    double coeff;
    double mainMeshSize;
    double otherMeshSize;
    GeometricFace* mainFace;
    GeometricFace* otherFace;
    // we use father as main mesh if it has smaller mesh size
    if (hFather <= hChild)
    {
        mainDomain = this;
        otherDomain = &child;
        coeff = 1;
        mainFace = &outlet;
        otherFace = &inlet;
        mainMeshSize = hFather;
        otherMeshSize = hChild;
        std::string msg = "Selecting father building block for interpolation\n";
        printlog(GREEN, msg, M_verbose);
    }
    else
    {
        mainDomain = &child;
        otherDomain = this;
        coeff = -1;
        mainFace = &inlet;
        otherFace = &outlet;
        mainMeshSize = hChild;
        otherMeshSize = hFather;
        std::string msg = "Selecting child building block for interpolation\n";
        printlog(GREEN, msg, M_verbose);
    }

    // create traces on both subdomains (to assemble the interpolation matrices)
    unsigned int mainNumTraces;
    VectorPtr* mainTraces;
    mainTraces = mainDomain->M_coupler.assembleTraces(*mainFace, 1.0,
                                                      mainDomain->M_couplingFESpace,
                                                      mainDomain->M_couplingFESpaceScalarETA,
                                                      mainNumTraces);

    unsigned int otherNumTraces;
    VectorPtr* otherTraces;
    otherTraces =
      otherDomain->M_coupler.assembleTraces(*otherFace, 1.0,
                                            otherDomain->M_couplingFESpace,
                                            otherDomain->M_couplingFESpaceScalarETA,
                                            otherNumTraces);

    otherDomain->M_coupler.buildInterpolators(M_datafile, mainMeshSize,
                                              otherMeshSize,
                                              mainDomain->getCouplingFESpace(),
                                              otherDomain->getCouplingFESpace(),
                                              mainFace, otherFace);

    otherDomain->M_coupler.buildInterpolationMatrices(mainTraces, mainNumTraces,
                                                      otherTraces, otherNumTraces);

    unsigned int nBasisFunctions = basisFunction->getNumBasisFunctions();

    VectorPtr* mainVectors;
    if (std::strcmp(basisFunction->getType().c_str(), "traces"))
        mainVectors = mainDomain->
                      M_coupler.assembleCouplingVectors(basisFunction,
                                                        *mainFace,
                                                        coeff,
                                                        mainDomain->M_couplingFESpace,
                                                        mainDomain->M_couplingFESpaceScalarETA,
                                                        nBasisFunctions);
    else
        mainVectors = mainDomain->
                      M_coupler.assembleTraces(*mainFace, coeff,
                                               mainDomain->M_couplingFESpace,
                                               mainDomain->M_couplingFESpaceScalarETA,
                                               nBasisFunctions, false);

    MapEpetraPtr lagrangeMultiplierMap = buildLagrangeMultiplierMap(nBasisFunctions,
                                                                    child,
                                                                    globalMap,
                                                                    maps,
                                                                    dimensions);

    MapEpetraPtr primalMap = mainDomain->M_primalMaps[mainDomain->M_indexCoupling];

    MatrixPtr couplingMatrix;
    couplingMatrix = mainDomain->M_coupler.fillMatrixWithVectorsInterpolated(mainVectors,
                                               nBasisFunctions,
                                               lagrangeMultiplierMap,
                                               primalMap,
                                               mainDomain->numberOfComponents(),
                                               otherDomain->M_treeNode->M_ID);

    MatrixPtr massMatrixMain = mainDomain->assembleBoundaryMatrix(*mainFace);
    MatrixPtr massMatrixOther = otherDomain->assembleBoundaryMatrix(*otherFace);

    mainDomain->M_coupler.buildCouplingMatrices(massMatrixMain,
                                                otherDomain->M_treeNode->M_ID);


    otherDomain->M_coupler.buildCouplingMatrices(massMatrixOther,
                                                 mainDomain->M_treeNode->M_ID,
                                                 couplingMatrix,
                                                 massMatrixMain);

    delete[] mainTraces;
    delete[] otherTraces;
}

AbstractAssembler::MapEpetraPtr
AbstractAssembler::
buildLagrangeMultiplierMap(const unsigned int nBasisFunctions,
                           AbstractAssembler& child,
                           AbstractAssembler::MapEpetraPtr& globalMap,
                           std::vector<AbstractAssembler::MapEpetraPtr>& maps,
                           std::vector<unsigned int>& dimensions)
{
    MapEpetraPtr lagrangeMultiplierMap;

    // build map for the lagrange multipliers (nBasisFunctions has been
    // modified in orthonormalization)
    unsigned int myel = (numberOfComponents() * nBasisFunctions) / M_comm->NumProc();
    // the first process takes care of the remainder
    if (M_comm->MyPID() == 0)
    {
        myel += (numberOfComponents() * nBasisFunctions) % M_comm->NumProc();
    }
    unsigned int mapSize = numberOfComponents() * nBasisFunctions;
    lagrangeMultiplierMap.reset(new
                    AbstractAssembler::MapEpetra(mapSize, myel,
                                                 0, M_comm));

    M_dualMaps.push_back(lagrangeMultiplierMap);
    child.M_dualMaps.push_back(lagrangeMultiplierMap);

    *globalMap += *lagrangeMultiplierMap;
    maps.push_back(lagrangeMultiplierMap);
    dimensions.push_back(lagrangeMultiplierMap->map(LifeV::Unique)
                                        ->NumGlobalElements());

    return lagrangeMultiplierMap;
}


void
AbstractAssembler::
assembleCouplingMatrices(AbstractAssembler& child,
                         const unsigned int& indexOutlet,
                         const unsigned int& interfaceIndex,
                         AbstractAssembler::MapEpetraPtr& globalMap,
                         std::vector<AbstractAssembler::MapEpetraPtr>& maps,
                         std::vector<unsigned int>& dimensions)
{
    std::string msg("[AbstractAssembler] start ");
    msg += "building coupling matrices between blocks with indices ";
    msg += std::to_string(M_treeNode->M_ID);
    msg += " and ";
    msg += std::to_string(child.M_treeNode->M_ID);
    msg += "\n";
    printlog(MAGENTA, msg, M_verbose);

    M_interfacesIndices.push_back(interfaceIndex);
    child.M_interfacesIndices.push_back(interfaceIndex);

    unsigned int nComponents = numberOfComponents();

    GeometricFace inlet = child.M_treeNode->M_block->getInlet();
    GeometricFace outlet = M_treeNode->M_block->getOutlet(indexOutlet);

    std::shared_ptr<BasisFunctionFunctor> basisFunction;
    basisFunction = castBasisFunction(inlet);
    unsigned int nBasisFunctions = basisFunction->getNumBasisFunctions();

    bool useInterpolation = M_datafile("coupling/use_interpolation", true);
    if (!std::strcmp(basisFunction->getType().c_str(), "traces") && !useInterpolation)
    {
        std::string msg;
        msg = "It is necessary to use interpolation ";
        msg +=  "when using trace basis functions!\n";
        throw Exception(msg);
    }

    if (useInterpolation)
    {
        assembleCouplingMatricesInterpolation(child, indexOutlet, interfaceIndex,
                                              basisFunction, globalMap,
                                              maps, dimensions);
        return;
    }

    VectorPtr* couplingVectorsFather;
    VectorPtr* couplingVectorsChild;
    couplingVectorsFather = M_coupler.assembleCouplingVectors(basisFunction,
                                                              outlet, 1,
                                                              M_couplingFESpace,
                                                              M_couplingFESpaceScalarETA,
                                                              nBasisFunctions);
    couplingVectorsChild = child.M_coupler.assembleCouplingVectors(basisFunction,
                                                                   inlet, -1,
                                                                   child.M_couplingFESpace,
                                                                   child.M_couplingFESpaceScalarETA,
                                                                   nBasisFunctions);

    MapEpetraPtr lagrangeMultiplierMap = buildLagrangeMultiplierMap(nBasisFunctions,
                                                                    child,
                                                                    globalMap,
                                                                    maps,
                                                                    dimensions);

    MatrixPtr massMatrixFather = assembleBoundaryMatrix(outlet);
    MatrixPtr massMatrixChild = child.assembleBoundaryMatrix(inlet);

    MapEpetraPtr primalMap = M_primalMaps[M_indexCoupling];

    M_coupler.fillMatrixWithVectorsInterpolated(couplingVectorsFather,
                                                nBasisFunctions,
                                                lagrangeMultiplierMap,
                                                primalMap,
                                                numberOfComponents(),
                                                child.M_treeNode->M_ID);

    primalMap = child.M_primalMaps[child.M_indexCoupling];

    child.M_coupler.fillMatrixWithVectorsInterpolated(couplingVectorsChild,
                                                      nBasisFunctions,
                                                      lagrangeMultiplierMap,
                                                      primalMap,
                                                      child.numberOfComponents(),
                                                      this->M_treeNode->M_ID);


    M_coupler.buildCouplingMatrices(massMatrixFather, child.M_treeNode->M_ID);
    child.M_coupler.buildCouplingMatrices(massMatrixChild, this->M_treeNode->M_ID);

    printlog(MAGENTA, "done\n", M_verbose);

    delete[] couplingVectorsFather;
    delete[] couplingVectorsChild;
}

AbstractAssembler::VectorPtr
AbstractAssembler::
reconstructLagrangeMultipliers(std::vector<VectorPtr> solutions, unsigned int offset)
{
    VectorPtr lagrangeMultipliers(new Vector(M_couplingFESpace->map()));
    lagrangeMultipliers->zero();
    unsigned int count = offset;
    for (std::map<unsigned int, MatrixPtr>::iterator it = M_coupler.getMapsQTsInterpolated().begin();
         it != M_coupler.getMapsQTsInterpolated().end(); it++)
    {
        // VectorPtr vcopy(new Vector(*solutions[count]));
        // vcopy->zero();
        // vcopy->operator[](7) = 1;
        // vcopy->operator[](12) = 1;
        // vcopy->operator[](18) = 1;
        // vcopy->operator[](23) = 1;
        // vcopy->operator[](30) = 1;
        // vcopy->operator[](33) = 1;
        // vcopy->operator[](37) = 1;
        // *lagrangeMultipliers += (*it->second) * (*vcopy);

        *lagrangeMultipliers += (*it->second) * (*solutions[count]);
        count++;
    }
    return lagrangeMultipliers;
}


// we copy the matrix in order to control boundary conditions
AbstractAssembler::MatrixPtr
AbstractAssembler::
getQT(const unsigned int& flag)
{
    MatrixPtr retMatrix(new Matrix(*M_coupler.getQT(flag)));
    // multiply by -1 because we are solving Hdu/dt = F!
    *retMatrix *= (-1);
    return retMatrix;
}

AbstractAssembler::MatrixPtr
AbstractAssembler::
getQ(const unsigned int& flag)
{
    MatrixPtr retMatrix(new Matrix(*M_coupler.getQ(flag)));
    // multiply by -1 because we are solving Hdu/dt = F!
    *retMatrix *= (-1);
    return retMatrix;
}

unsigned int
AbstractAssembler::
getIndexCoupling()
{
    return M_indexCoupling;
}

std::vector<unsigned int>
AbstractAssembler::
getInterfacesIndices()
{
    return M_interfacesIndices;
}

void
AbstractAssembler::
multiplyVectorsByMassMatrix(VectorPtr* couplingVectors,
                            const unsigned int& nBasisFunctions,
                            MatrixPtr massMatrix)
{
    for (int i = 0; i < nBasisFunctions; i++)
    {
        *couplingVectors[i] = (*massMatrix) * (*couplingVectors[i]);
    }
}

std::map<unsigned int, AbstractAssembler::MatrixPtr>&
AbstractAssembler::
getMapsQs()
{
    return M_coupler.getMapsQs();
}

std::map<unsigned int, AbstractAssembler::MatrixPtr>&
AbstractAssembler::
getMapsQTs()
{
    return M_coupler.getMapsQTs();
}

void
AbstractAssembler::
setTimestep(double dt)
{
    M_dt = dt;
}


}  // namespace RedMA