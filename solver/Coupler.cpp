#include <Coupler.hpp>

namespace RedMA
{

Coupler::
Coupler(commPtr_Type comm, bool verbose) :
  M_comm(comm),
  M_verbose(verbose)
{
}

void
Coupler::
gramSchmidt(Coupler::VectorPtr* basis1,
            Coupler::VectorPtr* basis2,
            Coupler::MatrixPtr massMatrix1,
            Coupler::MatrixPtr massMatrix2,
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
            double uv = Coupler::dotProd(basis1, basis2, i, j,
                                                   massMatrix1, massMatrix2);
            double uu = Coupler::dotProd(basis1, basis2, j, j,
                                                   massMatrix1, massMatrix2);

            *basis1[i] += (-uv/uu) * (*basis1[j]);
            *basis2[i] += (-uv/uu) * (*basis2[j]);
        }

        double curnorm = Coupler::dotProd(basis1, basis2, i, i,
                                                    massMatrix1, massMatrix2);
        *basis1[i] *= (1.0/std::sqrt(curnorm));
        *basis2[i] *= (1.0/std::sqrt(curnorm));
    }

    delete[] quickSearch;
    delete[] norms;
}

void
Coupler::
POD(VectorPtr*& basis1,
    VectorPtr*& basis2,
    MatrixPtr massMatrix1,
    MatrixPtr massMatrix2, unsigned int& nVectors,
    double tol,
    unsigned int offset)
{
    std::string msg = "[Coupler] performing POD ...\n";
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
            if (std::sqrt(eigenvalues[i]) < tol)
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

    msg = "[Coupler] done ...\n";
    printlog(MAGENTA, msg, M_verbose);
}

double
Coupler::
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

Coupler::VectorPtr*
Coupler::
assembleCouplingVectors(std::shared_ptr<BasisFunctionFunctor> basisFunction,
                        GeometricFace face, const double& coeff,
                        FESpacePtr couplingFespace,
                        ETFESpaceCouplingPtr couplingFESpaceETA,
                        unsigned int nBasisFunctions,
                        VectorPtr* otherInterfaceVectors,
                        InterpolationPtr interpolator)
{
    using namespace LifeV;
    // using namespace ExpressionAssembly;
    VectorPtr* couplingVectors = new VectorPtr[nBasisFunctions];
    MapEpetra couplingMap = couplingFESpaceETA->map();

    if (otherInterfaceVectors == nullptr && interpolator == nullptr)
    {
        unsigned int faceFlag = face.M_flag;
        MeshPtr mesh = couplingFESpaceETA->mesh();

        for (unsigned int i = 0; i < nBasisFunctions; i++)
        {
            // use repeated if you are integrating
            // VectorPtr currentMode(new Vector(couplingMap, Repeated));
            VectorPtr currentMode(new Vector(couplingMap, Unique));
            currentMode->zero();

            basisFunction->setIndex(i);

            // CoutRedirecter ct;
            // ct.redirect();

            BCFunctionBase bFunction(basisFunction->function());

            BoundaryConditionPtr bcs;
            bcs.reset(new BCHandler);

            bcs->addBC("Boundary_function", faceFlag, LifeV::Essential, LifeV::Full,
                        bFunction, 1);

            Function curFunction = basisFunction->function();
            couplingFespace->interpolateBC(*bcs, *currentMode, 0.0);
            *currentMode *= coeff;
            couplingVectors[i] = currentMode;
            // ct.restore();
        }
    }
    else
    {
        for (unsigned int i = 0; i < nBasisFunctions; i++)
        {
            VectorPtr currentMode(new Vector(couplingMap, Unique));
            currentMode->zero();

            VectorPtr otherVector3D(new Vector(couplingFespace->map()));
            // we need a 3D vector
            otherVector3D->subset(*otherInterfaceVectors[i],
                                 couplingMap,
                                 0,
                                 0);
            interpolator->updateRhs(otherVector3D);
            interpolator->interpolate();

            VectorPtr currentMode3D(new Vector(couplingFespace->map()));
            interpolator->solution(currentMode3D);
            currentMode->subset(*currentMode3D, couplingMap, 0, 0);

            *currentMode *= coeff;
            couplingVectors[i] = currentMode;
        }
    }
    return couplingVectors;
}

Coupler::VectorPtr*
Coupler::
assembleTraces(GeometricFace face, const double& coeff,
               FESpacePtr couplingFespace,
               ETFESpaceCouplingPtr couplingFESpaceETA,
               unsigned int& numBasisFunctions,
               bool keepZeros)
{
    // find indices of nodes on the interface
    unsigned int faceFlag = face.M_flag;
    MeshPtr mesh = couplingFESpaceETA->mesh();
    MapEpetra couplingMap = couplingFESpaceETA->map();

    VectorPtr indicatorFc(new Vector(couplingMap, LifeV::Unique));
    indicatorFc->zero();

    CoutRedirecter ct;
    ct.redirect();

    LifeV::BCFunctionBase bFunction(fOne);

    BoundaryConditionPtr bcs;
    bcs.reset(new LifeV::BCHandler);

    bcs->addBC("Boundary_function", faceFlag, LifeV::Essential, LifeV::Full,
                bFunction, 1);

    couplingFespace->interpolateBC(*bcs, *indicatorFc, 0.0);
    ct.restore();

    Epetra_Map primalMapEpetra = indicatorFc->epetraMap();
    unsigned int numElements = primalMapEpetra.NumMyElements();

    numBasisFunctions = 0;

    VectorPtr* couplingVectors;

    // add vectors also for the entries that are not on the boundary. Note
    // we are assembling scalar vectors even for 3D fields (because we are
    // using the map in FESpaceETA)
    if (keepZeros)
    {
        numBasisFunctions = numElements;
        couplingVectors = new VectorPtr[numBasisFunctions];

        for (unsigned int dof = 0; dof < numElements; dof++)
        {
            VectorPtr newVector(new Vector(*indicatorFc));
            newVector->zero();
            unsigned int gdof = primalMapEpetra.GID(dof);
            if (indicatorFc->isGlobalIDPresent(gdof))
            {
                double value(indicatorFc->operator[](gdof));
                if (std::abs(value-1) < 1e-15)
                {
                    newVector->operator[](gdof) = 1;
                }
            }
            couplingVectors[dof] = newVector;
            *couplingVectors[dof] *= coeff;
        }
    }
    else
    {
        for (unsigned int dof = 0; dof < numElements; dof++)
        {
            unsigned int gdof = primalMapEpetra.GID(dof);
            if (indicatorFc->isGlobalIDPresent(gdof))
            {
                double value(indicatorFc->operator[](gdof));
                if (std::abs(value-1) < 1e-15)
                {
                    numBasisFunctions++;
                }
            }
        }
        couplingVectors = new VectorPtr[numBasisFunctions];
        unsigned int count = 0;
        for (unsigned int dof = 0; dof < numElements; dof++)
        {
            unsigned int gdof = primalMapEpetra.GID(dof);
            if (indicatorFc->isGlobalIDPresent(gdof))
            {
                double value(indicatorFc->operator[](gdof));
                if (std::abs(value-1) < 1e-15)
                {
                    VectorPtr newVector(new Vector(*indicatorFc));
                    newVector->zero();
                    newVector->operator[](gdof) = 1;
                    couplingVectors[count] = newVector;
                    *couplingVectors[count] *= coeff;
                    count++;
                }
            }
        }
    }
    return couplingVectors;
}

Coupler::MatrixPtr
Coupler::
buildSingleInterpolationMatrix(VectorPtr* vectors,
                               const unsigned int numVectors,
                               MapEpetraPtr fromMap,
                               MapEpetraPtr toMap,
                               InterpolationPtr interpolator)
{
    using namespace LifeV;
    const Real dropTolerance(1e-14);

    MapEpetraPtr fromCompleteMap = interpolator->getKnownMap();
    MapEpetraPtr toCompleteMap = interpolator->getUnknownMap();

    VectorPtr* interpolatedVectors = new VectorPtr[numVectors];
    bool* indicesZeros = new bool[numVectors];
    for (unsigned int i = 0; i < numVectors; i++)
    {
        indicesZeros[i] = false;
        VectorPtr toVector(new Vector(*toMap, Unique));
        toVector->zero();

        VectorPtr fromVector3D(new Vector(*fromCompleteMap));
        // we need a 3D vector
        fromVector3D->subset(*vectors[i], *fromMap, 0, 0);
        if (vectors[i]->norm1() > 0)
        {
            interpolator->updateRhs(fromVector3D);
            interpolator->interpolate();

            VectorPtr toVector3D(new Vector(*toCompleteMap, Unique));
            interpolator->solution(toVector3D);
            toVector->subset(*toVector3D, *toMap, 0, 0);
        }
        else
        {
            indicesZeros[i] = true;
        }
        interpolatedVectors[i] = toVector;
    }

    MatrixPtr retMatrix(new Matrix(*toCompleteMap, numVectors, false));
    retMatrix->zero();

    Epetra_Map toMapEpetra = interpolatedVectors[0]->epetraMap();
    unsigned int toNumElements = toMapEpetra.NumMyElements();
    unsigned int toNDimensions = 3; // note: interpolation does not work for
                                    // non-threedimensional fields
    unsigned int toNTotalDofs = interpolatedVectors[0]->size();

    Epetra_Map fromMapEpetra = vectors[0]->epetraMap();
    unsigned int fromNumElements = fromMapEpetra.NumMyElements();
    unsigned int fromNDimensions = 3;
    unsigned int fromNTotalDofs = vectors[0]->size();

    for (unsigned int dim = 0; dim < toNDimensions; dim++)
    {
        for (unsigned int i = 0; i < numVectors; i++)
        {
            if (!indicesZeros[i])
            {
                Vector vectorUnique(*interpolatedVectors[i], Unique);
                for (unsigned int dof = 0; dof < toNumElements; dof++)
                {
                    unsigned int gdof = toMapEpetra.GID(dof);
                    if (vectorUnique.isGlobalIDPresent(gdof))
                    {
                        double value(vectorUnique[gdof]);
                        if (std::abs(value) > dropTolerance)
                        {
                            retMatrix->addToCoefficient(gdof + dim * toNTotalDofs,
                                                        i + dim * fromNTotalDofs,
                                                        value);
                        }
                    }
                }
            }
        }
    }

    M_comm->Barrier();
    retMatrix->globalAssemble(fromCompleteMap, toCompleteMap);

    delete[] indicesZeros;
    return retMatrix;
}


void
Coupler::
buildInterpolationMatrices(VectorPtr* mainTraces, const unsigned int mainNumTraces,
                           VectorPtr* otherTraces, const unsigned int otherNumTraces)
{
    MapEpetraPtr mainMap =  mainTraces[0]->mapPtr();
    MapEpetraPtr otherMap =  otherTraces[0]->mapPtr();

    M_matrixInterpolationMainToOther =
                      buildSingleInterpolationMatrix(mainTraces, mainNumTraces,
                                                     mainMap, otherMap,
                                                     M_interpolatorMainToOther);
    M_matrixInterpolationOtherToMain =
                      buildSingleInterpolationMatrix(otherTraces, otherNumTraces,
                                                     otherMap, mainMap,
                                                     M_interpolatorOtherToMain);
}

void
Coupler::
buildInterpolators(GetPot datafile, double mainMeshSize, double otherMeshSize,
                   FESpacePtr mainFESpace, FESpacePtr otherFESpace,
                   GeometricFace* mainFace, GeometricFace* otherFace)
{
    M_interpolatorMainToOther = buildSingleInterpolator(datafile, mainMeshSize,
                                                 mainFESpace, otherFESpace,
                                                 mainFace, otherFace);
    M_interpolatorOtherToMain = buildSingleInterpolator(datafile, otherMeshSize,
                                                 otherFESpace, mainFESpace,
                                                 otherFace, mainFace);
}

Coupler::InterpolationPtr
Coupler::
buildSingleInterpolator(GetPot datafile, double meshSize,
                        FESpacePtr FESpaceFrom,
                        FESpacePtr FESpaceTo,
                        GeometricFace* faceFrom,
                        GeometricFace* faceTo)
{
    Teuchos::RCP<Teuchos::ParameterList> belosList =
                                  Teuchos::rcp (new Teuchos::ParameterList);
    belosList = Teuchos::getParametersFromXmlFile("SolverParamList_rbf3d.xml");

    InterpolationPtr interpolator;
    interpolator.reset(new Interpolation);
    interpolator->setup(datafile, belosList);
    interpolator->setMeshSize(meshSize);

    VectorPtr dummyVectorFrom(new
                           Vector(FESpaceFrom->map()));
    dummyVectorFrom->zero();

    VectorPtr dummyVectorTo(new
                          Vector(FESpaceTo->map()));
    dummyVectorTo->zero();

    interpolator->setVectors(dummyVectorFrom, dummyVectorTo);
    interpolator->setFlag(faceFrom->M_flag);
    interpolator->buildTableDofs_known(FESpaceFrom);
    interpolator->identifyNodes_known();
    interpolator->buildKnownInterfaceMap();

    interpolator->setFlag(faceTo->M_flag);
    interpolator->buildTableDofs_unknown(FESpaceTo);
    interpolator->identifyNodes_unknown();
    interpolator->buildUnknownInterfaceMap();

    interpolator->setFlag(faceFrom->M_flag);
    interpolator->buildOperators();

    return interpolator;
}

double
Coupler::
fOne(const double& t, const double& x, const double& y, const double& z,
     const unsigned int& i)
{
    return 1.0;
}


void
Coupler::
fillMatricesWithVectors(VectorPtr* couplingVectors,
                        const unsigned int& nBasisFunctions,
                        MapEpetraPtr lagrangeMap,
                        MapEpetraPtr map,
                        unsigned int numberOfComponents,
                        const unsigned int& flagAdjacentDomain)
{
    using namespace LifeV;
    const Real dropTolerance(2.0 * std::numeric_limits<Real>::min());
    unsigned faceFlag = flagAdjacentDomain;
    // note:: we have to specify the second argument of the constructor (number
    // of elements per row)
    MatrixPtr QT(new Matrix(*map, nBasisFunctions, false));
    QT->zero();

    MatrixPtr Q(new Matrix(*lagrangeMap));
    Q->zero();

    Epetra_Map primalMapEpetra = couplingVectors[0]->epetraMap();
    unsigned int numElements = primalMapEpetra.NumMyElements();
    // unsigned int nTotalDofs = primalFespace->dof().numTotalDof();
    unsigned int nTotalDofs = couplingVectors[0]->size();
    for (unsigned int dim = 0; dim < numberOfComponents; dim++)
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
    M_comm->Barrier();
    QT->globalAssemble(lagrangeMap, map);

    Q = QT->transpose();
    Q->globalAssemble(map, lagrangeMap);

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
}

Coupler::MatrixPtr
Coupler::
fillMatrixWithVectorsInterpolated(VectorPtr* couplingVectors,
                                  const unsigned int& nBasisFunctions,
                                  MapEpetraPtr lagrangeMap,
                                  MapEpetraPtr map,
                                  unsigned int numberOfComponents,
                                  const unsigned int& flagAdjacentDomain)
{
    using namespace LifeV;
    const Real dropTolerance(2.0 * std::numeric_limits<Real>::min());
    unsigned faceFlag = flagAdjacentDomain;
    // note:: we have to specify the second argument of the constructor (number
    // of elements per row)
    MatrixPtr QT(new Matrix(*map, nBasisFunctions, false));
    QT->zero();

    Epetra_Map primalMapEpetra = couplingVectors[0]->epetraMap();
    unsigned int numElements = primalMapEpetra.NumMyElements();
    // unsigned int nTotalDofs = primalFespace->dof().numTotalDof();
    unsigned int nTotalDofs = couplingVectors[0]->size();
    for (unsigned int dim = 0; dim < numberOfComponents; dim++)
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
    M_comm->Barrier();
    QT->globalAssemble(lagrangeMap, map);

    if (M_mapQTsInterpolated.find(faceFlag) == M_mapQTsInterpolated.end())
    {
        M_mapQTsInterpolated[faceFlag] = QT;
    }
    else
    {
        std::string errorMsg("Coupling matrix interpolated with key = ");
        errorMsg += std::to_string(faceFlag) + " have already been assembled!\n";

        throw Exception(errorMsg);
    }
    return QT;
}

void
Coupler::
buildCouplingMatrices(MatrixPtr myMass,
                      const unsigned int& flagAdjacentDomain,
                      MatrixPtr matrixToInterpolate,
                      MatrixPtr otherMass)
{
    MatrixPtr Qs;
    MatrixPtr QTs;
    // then we should have already computed our own coupling matrices and we
    // can just multiply by the mass matrix
    if (matrixToInterpolate == nullptr)
    {
        MatrixPtr QTsInt = M_mapQTsInterpolated[flagAdjacentDomain];

        MapEpetra map1 = QTsInt->rangeMap();
        MapEpetra map2 = QTsInt->domainMap();

        Qs.reset(new Matrix(map2));
        QTs.reset(new Matrix(map1));

        myMass->multiply(false, *QTsInt, false, *QTs, true);
        QTsInt->multiply(true, *myMass, false, *Qs, true);
    }
    // we need to use the interpolation matrices. We suppose that we are on the
    // "other" domain (the one towards which we interpolate)
    else
    {
        MapEpetra map = myMass->rangeMap();

        // first we handle the transpose matrix
        QTs.reset(new Matrix(map));
        MatrixPtr res(new Matrix(map));

        M_matrixInterpolationMainToOther->multiply(false, *matrixToInterpolate,
                                                   false, *res, true);

        myMass->multiply(false, *res, false, *QTs, true);

        // here we put the res as our own matrix (for exporting solution)
        M_mapQTsInterpolated[flagAdjacentDomain] =  res;
        *M_mapQTsInterpolated[flagAdjacentDomain] *= (-1.0);

        // handle non transpose matrix
        map = otherMass->rangeMap();
        res.reset(new Matrix(map));

        otherMass->multiply(false, *M_matrixInterpolationOtherToMain, false,
                            *res, true);
        // otherMass->multiply(false, *M_matrixInterpolationMainToOther, true,
        //                     *res, true);

        map = matrixToInterpolate->domainMap();
        Qs.reset(new Matrix(map));

        matrixToInterpolate->multiply(true, *res, false, *Qs, true);
        *Qs *= (-1.0);
        *QTs *= (-1.0);
    }
    M_mapQTs[flagAdjacentDomain] = QTs;
    M_mapQs[flagAdjacentDomain] = Qs;
}

Coupler::MatrixPtr
Coupler::
getQT(const unsigned int& flag)
{
    return M_mapQTs[flag];
}

Coupler::MatrixPtr
Coupler::
getQ(const unsigned int& flag)
{
    return M_mapQs[flag];
}

double*
Coupler::
computeCorrelationMatrix(Coupler::VectorPtr* basis1,
                         Coupler::VectorPtr* basis2,
                         Coupler::MatrixPtr mass1,
                         Coupler::MatrixPtr mass2,
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

std::map<unsigned int, Coupler::MatrixPtr>&
Coupler::
getMapsQTs()
{
    return M_mapQTs;
}

std::map<unsigned int, Coupler::MatrixPtr>&
Coupler::
getMapsQTsInterpolated()
{
    return M_mapQTsInterpolated;
}

std::map<unsigned int, Coupler::MatrixPtr>&
Coupler::
getMapsQs()
{
    return M_mapQs;
}

}  // namespace RedMA
