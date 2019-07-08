#include <AbstractAssembler.hpp>

namespace RedMA
{

AbstractAssembler::
AbstractAssembler(const GetPot& datafile, commPtr_Type comm,
                  const TreeNodePtr& treeNode, bool verbose) :
  M_datafile(datafile),
  M_comm(comm),
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
gramSchmidt(AbstractAssembler::VectorPtr* basis1,
            AbstractAssembler::VectorPtr* basis2,
            AbstractAssembler::MatrixPtr massMatrix1,
            AbstractAssembler::MatrixPtr massMatrix2,
            unsigned int& nVectors)
{
    for (unsigned int i = 0; i < nVectors; i++)
    {
        VectorPtr v1(new Vector(*basis1[i], LifeV::Unique));
        VectorPtr v2(new Vector(*basis2[i], LifeV::Unique));
        basis1[i] = v1;
        basis2[i] = v2;
    }

    bool* quickSearch = new bool[nVectors];
    unsigned int vectorsZeroed = 0;
    double* norms = new double[nVectors];
    for (unsigned int i = 0; i < nVectors; i++)
    {
        double u1 = std::sqrt(dotProd(basis1, basis2, i, i));
        norms[i] = u1;
        quickSearch[i] = false;
        if (u1 < 5e-12)
        {
            quickSearch[i] = true;
            basis1[i]->zero();
            basis2[i]->zero();
            vectorsZeroed++;
        }
        else
        {
            for (unsigned int j = 0; j < i; j++)
            {
                if (!quickSearch[j])
                {
                    double u2 = dotProd(basis1, basis2, i, j);
                    double ratio = std::abs(u2)/(norms[i]*norms[j]);
                    double theta = std::acos(ratio);
                    if (std::abs(theta) < 5e-1)
                    {
                        quickSearch[i] = true;
                        basis1[i]->zero();
                        basis2[i]->zero();
                        vectorsZeroed++;
                        break;
                    }
                }
            }
        }
    }
    unsigned int auxIndex = 0;

    for (unsigned int i = 0; i < nVectors - vectorsZeroed; i++)
    {
        while (quickSearch[auxIndex] && auxIndex < nVectors)
        {
            auxIndex++;
        }
        basis1[i] = basis1[auxIndex];
        basis2[i] = basis2[auxIndex];
        auxIndex++;
    }

    nVectors -= vectorsZeroed;

    for (unsigned int i = 0; i < nVectors; i++)
    {
        for (unsigned int j = 0; j < i; j++)
        {
            double uv = AbstractAssembler::dotProd(basis1, basis2, i, j,
                                                   massMatrix1, massMatrix2);
            double uu = AbstractAssembler::dotProd(basis1, basis2, j, j,
                                                   massMatrix1, massMatrix2);

            *basis1[i] += (-uv/uu) * (*basis1[j]);
            *basis2[i] += (-uv/uu) * (*basis2[j]);
        }

        double curnorm = AbstractAssembler::dotProd(basis1, basis2, i, i,
                                                    massMatrix1, massMatrix2);
        *basis1[i] *= (1.0/std::sqrt(curnorm));
        *basis2[i] *= (1.0/std::sqrt(curnorm));
    }

    delete[] quickSearch;
    delete[] norms;
}

void
AbstractAssembler::
POD(VectorPtr*& basis1,
    VectorPtr*& basis2,
    MatrixPtr massMatrix1,
    MatrixPtr massMatrix2, unsigned int& nVectors,
    double tol,
    unsigned int offset)
{
    std::string msg = "[AbstractAssembler] performing POD ...\n";
    printlog(MAGENTA, msg, M_verbose);
    // subtracting constant
    nVectors = nVectors - offset;
    double* correlationMatrix =
            computeCorrelationMatrix(basis1 + offset, basis2 + offset, massMatrix1, massMatrix2,
                                     nVectors);

    double* fieldVL;
    double* eigenvalues;

    unsigned int initialN = nVectors;
    unsigned int length = nVectors * nVectors;

    if (M_comm->MyPID() == 0)
    {
        Epetra_LAPACK lapack;

        // double* singularValues = new double[M_nSnapshots];

        eigenvalues = new double[nVectors];
        fieldVL = new double[length];
        double* fieldVR = new double[1];

        int info;

        int lwork(-1);
        double* work(nullptr);
        work = new double[1];

        fieldVL = new double[length];
        fieldVR = new double[1];

        // compute optimal work in lwork
        lapack.GESVD('A', 'N', nVectors, nVectors, correlationMatrix, nVectors,
                      eigenvalues, fieldVL, nVectors, fieldVR, 1, work, &lwork,
                      &info);

        lwork = work[0];
        work = new double[lwork];

        lapack.GESVD('A', 'N', nVectors, nVectors, correlationMatrix, nVectors,
                      eigenvalues, fieldVL, nVectors, fieldVR, 1, work, &lwork,
                      &info);

        unsigned int i;
        // compute new number of nVectors
        for (i = 0; i < nVectors; i++)
        {
            std::string eigMsg = std::to_string(i) + ": " +
                                 std::to_string(std::sqrt(eigenvalues[i])) + "\n";
            printlog(GREEN, eigMsg, M_verbose);
            if (eigenvalues[i] < tol)
                break;
        }

        nVectors = i;
    }

    int nVectorsInt = nVectors;
    M_comm->Broadcast(&nVectorsInt, 1, 0);
    M_comm->Barrier();
    nVectors = nVectorsInt;

    if (M_comm->MyPID() != 0)
    {
        fieldVL = new double[length];
        eigenvalues = new double[nVectors];
    }

    M_comm->Broadcast(eigenvalues, nVectors, 0);
    M_comm->Broadcast(fieldVL, length, 0);
    M_comm->Barrier();

    // re-adding constant
    VectorPtr* newVectors1 = new VectorPtr[nVectors + offset];
    VectorPtr* newVectors2 = new VectorPtr[nVectors + offset];

    for (unsigned int i = 0; i < offset; i++)
    {
        newVectors1[i].reset(new Vector(basis1[0]->map(), LifeV::Unique));
        newVectors2[i].reset(new Vector(basis2[0]->map(), LifeV::Unique));
        *newVectors1[i] = *basis1[i];
        *newVectors2[i] = *basis2[i];
    }

    for (unsigned int i = 0; i < nVectors; i++)
    {
        double coeff = 1.0 / std::sqrt(eigenvalues[i]);

        newVectors1[i+offset].reset(new Vector(basis1[0]->map(), LifeV::Unique));
        newVectors2[i+offset].reset(new Vector(basis2[0]->map(), LifeV::Unique));
        newVectors1[i+offset]->zero();
        newVectors2[i+offset]->zero();

        for (unsigned int j = 0; j < initialN; j++)
        {
            unsigned int index = i * initialN + j;

            *newVectors1[i+offset] += (coeff * fieldVL[index]) * (*basis1[j+offset]);
            *newVectors2[i+offset] += (coeff * fieldVL[index]) * (*basis2[j+offset]);
        }
    }
    delete[] basis1;
    delete[] basis2;

    basis1 = newVectors1;
    basis2 = newVectors2;
    nVectors = nVectors + offset;

    msg = "[AbstractAssembler] done ...\n";
    printlog(MAGENTA, msg, M_verbose);
}

double
AbstractAssembler::
dotProd(VectorPtr* basis1, VectorPtr* basis2, unsigned int index1,
        unsigned int index2, MatrixPtr mass1, MatrixPtr mass2)
{
    double prod1 = 0, prod2 = 0;
    if (mass1 && mass2)
    {
        VectorPtr aux1(new Vector(basis1[index1]->map()));
        VectorPtr aux2(new Vector(basis2[index1]->map()));
        *aux1 = (*mass1) * (*basis1[index1]);
        *aux2 = (*mass2) * (*basis2[index1]);
        aux1->dot(*basis1[index2], prod1);
        aux2->dot(*basis2[index2], prod2);
    }
    else
    {
        basis1[index1]->dot(*basis1[index2], prod1);
        basis2[index1]->dot(*basis2[index2], prod2);
    }
    return prod1 + prod2;
}

AbstractAssembler::VectorPtr*
AbstractAssembler::
assembleCouplingVectors(std::shared_ptr<BasisFunctionFunctor> basisFunction,
                        GeometricFace face, const double& coeff)
{
    using namespace LifeV;
    // using namespace ExpressionAssembly;

    QuadratureBoundary boundaryQuadRule(buildTetraBDQR(quadRuleTria7pt));

    unsigned int nBasisFunctions = basisFunction->getNumBasisFunctions();

    VectorPtr* couplingVectors = new VectorPtr[nBasisFunctions];
    MapEpetra couplingMap = M_couplingFESpaceETA->map();

    unsigned int faceFlag = face.M_flag;
    MeshPtr mesh = M_couplingFESpaceETA->mesh();

    for (unsigned int i = 0; i < nBasisFunctions; i++)
    {
        // use repeated if you are integrated
        // VectorPtr currentMode(new Vector(couplingMap, Repeated));
        VectorPtr currentMode(new Vector(couplingMap, Unique));
        basisFunction->setIndex(i);

        CoutRedirecter ct;
        ct.redirect();

        BCFunctionBase bFunction(basisFunction->function());

        BoundaryConditionPtr bcs;
        bcs.reset(new BCHandler);

        bcs->addBC("Boundary_function", faceFlag, LifeV::Essential, LifeV::Full,
                    bFunction, 1);

        Function curFunction = basisFunction->function();
        M_couplingFESpace->interpolateBC(*bcs, *currentMode, 0.0);
        *currentMode *= coeff;
        couplingVectors[i] = currentMode;
        ct.restore();
        // integrate(boundary(mesh, faceFlag),
        //           boundaryQuadRule,
        //           M_couplingFESpaceETA,
        //           value(coeff) * eval(basisFunction, X) * phi_i
        //       ) >> currentMode;
        // couplingVectors[i] = currentMode;
    }
    M_couplingVector.reset(new Vector(*couplingVectors[1]));
    return couplingVectors;
}


AbstractAssembler::MatrixPtr
AbstractAssembler::
assembleBoundaryMatrix(GeometricFace face)
{
    using namespace LifeV;
    using namespace ExpressionAssembly;

    QuadratureBoundary boundaryQuadRule(buildTetraBDQR(quadRuleTria7pt));

    unsigned int faceFlag = face.M_flag;
    MeshPtr mesh = M_couplingFESpaceETA->mesh();

    MapEpetra couplingMap = M_couplingFESpaceETA->map();
    MatrixPtr boundaryMassMatrix(new Matrix(couplingMap));

    // assemble boundary mass matrix to orthonormalize w.r.t. L2 product
    integrate(boundary(mesh, faceFlag),
              boundaryQuadRule,
              M_couplingFESpaceETA,
              M_couplingFESpaceETA,
              phi_i * phi_j
              ) >> boundaryMassMatrix;
    boundaryMassMatrix->globalAssemble();

    return boundaryMassMatrix;
}

void
AbstractAssembler::
fillMatricesWithVectors(VectorPtr* couplingVectors,
                        const unsigned int& nBasisFunctions,
                        MapEpetraPtr lagrangeMap,
                        const unsigned int& flagAdjacentDomain)
{
    using namespace LifeV;
    const Real dropTolerance(2.0 * std::numeric_limits<Real>::min());
    unsigned faceFlag = flagAdjacentDomain;
    // note:: we have to specify the second argument of the constructor (number
    // of elements per row)
    MatrixPtr QT(new Matrix(*M_primalMaps[M_indexCoupling],
                            nBasisFunctions, false));
    QT->zero();

    MatrixPtr Q(new Matrix(*lagrangeMap));
    Q->zero();

    Epetra_Map primalMapEpetra = couplingVectors[0]->epetraMap();
    unsigned int numElements = primalMapEpetra.NumMyElements();
    // unsigned int nTotalDofs = primalFespace->dof().numTotalDof();
    unsigned int nTotalDofs = couplingVectors[0]->size();
    for (unsigned int dim = 0; dim < numberOfComponents(); dim++)
    {
        for (unsigned int i = 0; i < nBasisFunctions; i++)
        {
            Vector couplingVectorUnique(*couplingVectors[i], Unique);
            for (unsigned int dof = 0; dof < numElements; dof++)
            {
                unsigned int gdof = primalMapEpetra.GID(dof);
                if (couplingVectorUnique.isGlobalIDPresent(gdof))
                {
                    double value(couplingVectorUnique[gdof]);
                    if (std::abs(value) > dropTolerance)
                    {
                        // this must be changed to insert coefficients
                        QT->addToCoefficient(gdof + dim * nTotalDofs,
                                             i + dim * nBasisFunctions,
                                             value);
                    }
                }
            }
        }
    }
    QT->globalAssemble(lagrangeMap, M_primalMaps[M_indexCoupling]);

    Q = QT->transpose();
    Q->globalAssemble(M_primalMaps[M_indexCoupling], lagrangeMap);

    if (M_mapQTs.find(faceFlag) == M_mapQTs.end() &&
        M_mapQs.find(faceFlag) == M_mapQs.end())
    {
        M_mapQTs[faceFlag] = QT;
        M_mapQs[faceFlag] = Q;
    }
    else
    {
        std::string errorMsg("Coupling matrices with key = ");
        errorMsg += std::to_string(faceFlag) + " have already been assembled!\n";

        throw Exception(errorMsg);
    }

    delete[] couplingVectors;
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

    unsigned int nComponents = numberOfComponents();
    unsigned int nBasisFunctions;
    MapEpetraPtr lagrangeMultiplierMap;

    std::string typeBasis = M_datafile("coupling/type", "fourier");


    std::shared_ptr<BasisFunctionFunctor> basisFunction;

    GeometricFace inlet = child.M_treeNode->M_block->getInlet();
    GeometricFace outlet = M_treeNode->M_block->getOutlet(indexOutlet);

    if (!std::strcmp(typeBasis.c_str(), "fourier"))
    {
        unsigned int frequenciesTheta = M_datafile("coupling/frequencies_theta", 1);
        unsigned int frequenciesRadial = M_datafile("coupling/frequencies_radial", 1);

        basisFunction.reset(new FourierBasisFunction(inlet, frequenciesTheta,
                                                     frequenciesRadial));
    }
    else if (!std::strcmp(typeBasis.c_str(), "zernike"))
    {
        unsigned int nMax = M_datafile("coupling/nMax", 5);
        basisFunction.reset(new ZernikeBasisFunction(inlet, nMax));
    }

    nBasisFunctions = basisFunction->getNumBasisFunctions();
    VectorPtr* couplingVectorsFather =
                assembleCouplingVectors(basisFunction, outlet, 1);
    MatrixPtr massMatrixFather = assembleBoundaryMatrix(outlet);

    VectorPtr* couplingVectorsChild =
          child.assembleCouplingVectors(basisFunction, inlet, -1);
    MatrixPtr massMatrixChild = child.assembleBoundaryMatrix(inlet);
    unsigned int prev = nBasisFunctions;

    std::string orthStrategy = M_datafile("coupling/orthonormalization", "POD");

    if (std::strcmp(orthStrategy.c_str(), "none"))
    {
        if (!std::strcmp(orthStrategy.c_str(), "POD"))
        {
            double tol = M_datafile("coupling/beta_threshold", 9e-1);

            POD(couplingVectorsFather, couplingVectorsChild, massMatrixFather,
                massMatrixChild, nBasisFunctions, tol, 1);
            POD(couplingVectorsFather, couplingVectorsChild, massMatrixFather,
                massMatrixChild, nBasisFunctions, 0);
        }
        else if (!std::strcmp(orthStrategy.c_str(), "gram_schmidt"))
            gramSchmidt(couplingVectorsFather, couplingVectorsChild, massMatrixFather,
                        massMatrixChild, nBasisFunctions);
        else
        {
            std::string errMsg = "Orthonormalization method " + orthStrategy;
            errMsg += " does not exist!\n";

            throw Exception(errMsg);
        }
        msg = "Orthonormalizing with " + orthStrategy;
        msg += ": reducing number of bfs from " + std::to_string(prev) + " to " +
               std::to_string(nBasisFunctions) + "\n";
        printlog(GREEN, msg, M_verbose);
    }

    // up to this moment the coupling vectors contain only evaluation of basis functions
    // at the nodes. We need to multiply by the boundary mass matrices
    multiplyVectorsByMassMatrix(couplingVectorsFather, nBasisFunctions,
                                massMatrixFather);

    multiplyVectorsByMassMatrix(couplingVectorsChild, nBasisFunctions,
                                massMatrixChild);

    // build map for the lagrange multipliers (nBasisFunctions has been
    // modified in orthonormalization)
    unsigned int myel = (nComponents * nBasisFunctions) / M_comm->NumProc();
    // the first process takes care of the remainder
    if (M_comm->MyPID() == 0)
    {
        myel += (nComponents * nBasisFunctions) % M_comm->NumProc();
    }
    unsigned int mapSize = nComponents * nBasisFunctions;
    lagrangeMultiplierMap.reset(new
                    AbstractAssembler::MapEpetra(mapSize, myel,
                                                 0, M_comm));
    M_dualMaps.push_back(lagrangeMultiplierMap);
    child.M_dualMaps.push_back(lagrangeMultiplierMap);

    fillMatricesWithVectors(couplingVectorsFather, nBasisFunctions,
                            lagrangeMultiplierMap,
                            child.M_treeNode->M_ID);

    child.fillMatricesWithVectors(couplingVectorsChild, nBasisFunctions,
                                  lagrangeMultiplierMap, M_treeNode->M_ID);

    *globalMap += *lagrangeMultiplierMap;
    maps.push_back(lagrangeMultiplierMap);
    dimensions.push_back(lagrangeMultiplierMap->map(LifeV::Unique)
                                              ->NumGlobalElements());
    M_interfacesIndices.push_back(interfaceIndex);
    child.M_interfacesIndices.push_back(interfaceIndex);
    printlog(MAGENTA, "done\n", M_verbose);
}

// we copy the matrix in order to control boundary conditions
AbstractAssembler::MatrixPtr
AbstractAssembler::
getQT(const unsigned int& flag)
{
    MatrixPtr retMatrix(new Matrix(*M_mapQTs[flag]));
    // multiply by -1 because we are solving Hdu/dt = F!
    *retMatrix *= (-1);
    return retMatrix;
}

AbstractAssembler::MatrixPtr
AbstractAssembler::
getQ(const unsigned int& flag)
{
    MatrixPtr retMatrix(new Matrix(*M_mapQs[flag]));
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

double*
AbstractAssembler::
computeCorrelationMatrix(AbstractAssembler::VectorPtr* basis1,
                         AbstractAssembler::VectorPtr* basis2,
                         AbstractAssembler::MatrixPtr mass1,
                         AbstractAssembler::MatrixPtr mass2,
                         const unsigned int& nBasisFunctions)
{
    unsigned int matSize = nBasisFunctions * nBasisFunctions;

    for (unsigned int i = 0; i < nBasisFunctions; i++)
    {
        VectorPtr v1(new Vector(*basis1[i], LifeV::Unique));
        VectorPtr v2(new Vector(*basis2[i], LifeV::Unique));
        basis1[i] = v1;
        basis2[i] = v2;
    }

    double* correlationMatrix = new double[matSize];

    for (unsigned int i = 0; i < nBasisFunctions; i++)
    {
        for (unsigned int j = 0; j <= i; j++)
        {
            double prod = dotProd(basis1, basis2, i, j, mass1, mass2);
            correlationMatrix[nBasisFunctions * i + j] = prod;
            correlationMatrix[nBasisFunctions * j + i] = prod;
        }
    }

    return correlationMatrix;
}

void
AbstractAssembler::
multiplyVectorsByMassMatrix(VectorPtr* couplingVectors,
                            const unsigned int& nBasisFunctions,
                            MatrixPtr massMatrix)
{
    // std::cout << "------------" << std::endl;
    for (int i = 0; i < nBasisFunctions; i++)
    {
        Vector c(*couplingVectors[i]);
        *couplingVectors[i] = (*massMatrix) * (*couplingVectors[i]);
        // double prod;
        // c.dot(*couplingVectors[i], prod);
        // std::cout << "prod = " << prod << std::endl;
    }
}

}  // namespace RedMA
