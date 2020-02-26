#include "MDEIM.hpp"

namespace RedMA
{

MDEIM::
MDEIM()
{

}

void
MDEIM::
setDataContainer(const DataContainer& dataContainer)
{
    M_data = dataContainer;
}

void
MDEIM::
setComm(EPETRACOMM comm)
{
    M_comm = comm;
}

void
MDEIM::
setAssembler(SHP(aAssembler<FEVECTOR COMMA FEMATRIX>) assembler)
{
    M_assembler = assembler;
}

void
MDEIM::
addSnapshot(BlockMatrix<MatrixEp> newSnapshot)
{
    if (M_snapshots.size() == 0)
    {
        initializeMDEIMStructures(newSnapshot);
        M_nRows = newSnapshot.nRows();
        M_nCols = newSnapshot.nCols();
    }

    M_snapshots.push_back(newSnapshot);
}

void
MDEIM::
initializeMDEIMStructures(BlockMatrix<MatrixEp> matrix)
{
    unsigned int nrows = matrix.nRows();
    unsigned int ncols = matrix.nCols();
    M_structures.resize(nrows, ncols);

    for (unsigned int i = 0; i < nrows; i++)
    {
        for (unsigned int j = 0; j < ncols; j++)
        {
            initializeSingleMDEIMStructure(i, j, matrix.block(i, j));
        }
    }
}

void
MDEIM::
initializeSingleMDEIMStructure(const unsigned int& i,
                               const unsigned int& j, MatrixEp matrix)
{
    SHP(MATRIXEPETRA) matrixEpetra = matrix.data();

    if (matrixEpetra)
    {
        SHP(SingleMDEIMStructure) ms(new SingleMDEIMStructure());

        ms->numGlobalNonzeros = matrixEpetra->matrixPtr()->NumGlobalNonzeros();
        ms->numMyNonzeros = matrixEpetra->matrixPtr()->NumMyNonzeros();
        ms->numMyRows = matrixEpetra->matrixPtr()->NumMyRows();

        ms->numMyEntries = new int[ms->numMyRows];
        ms->columnIndeces = new int* [ms->numMyRows];
        int numEntries;

        for (unsigned int iR = 0; iR < ms->numMyRows; iR++)
        {
            ms->numMyEntries[iR] = matrixEpetra->matrixPtr()->NumMyEntries(iR);
            ms->columnIndeces[iR] = new int[ms->numMyEntries[iR]];
            double * values = new double[ms->numMyEntries[iR]];
            matrixEpetra->matrixPtr()->ExtractMyRowCopy(iR, ms->numMyEntries[iR],
                                                        numEntries, values,
                                                        ms->columnIndeces[iR]);

            for(unsigned int iC; iC < numEntries; iC++)
            {
                ms->columnIndeces[iR][iC] = matrixEpetra->matrixPtr()->ColMap().
                                                    GID(ms->columnIndeces[iR][iC]);
            }
            delete [] values;
        }

        ms->partialSumMyEntries = new int[ms->numMyRows + 1];
        ms->partialSumMyEntries[0] = 0;

        for (unsigned int iR = 0; iR < ms->numMyRows; iR++)
        {
            ms->partialSumMyEntries[iR+1] = ms->partialSumMyEntries[iR] + ms->numMyEntries[iR];
        }

        int * myGlobalElements = 0;

        ms->vectorMap.reset(new MAPEPETRA(ms->numGlobalNonzeros,
                                         ms->numMyNonzeros,
                                         myGlobalElements, M_comm));

        ms->allocated = true;
        M_structures(i,j) = ms;
    }
}

void
MDEIM::
performMDEIM()
{
    vectorizeSnapshots();
    performPOD();
}

void
MDEIM::
pickMagicPoints(const unsigned int& i, const unsigned int& j)
{
    std::vector<SHP(VECTOREPETRA)> basis = M_bases(i,j);
    SHP(SingleMDEIMStructure) ms = M_structures(i,j);
    // get first coefficient manually
    int maxInfo[3];
    Epetra_SerialDenseVector interpCoef(1);

    unsigned int N = basis.size();

    ms->N = N;
    ms->Qj.Reshape(1,1);

    // I keep these structures even if I don't plan to call this routine in parallel
    ms->localIndecesMagicPoints.resize(N);
    ms->magicPointsProcOwner.resize(N);
    ms->globalIndecesMagicPoints.resize(N);

    rbLifeV::getInfNorm(*(basis[0]), maxInfo);

    ms->localIndecesMagicPoints[0] = maxInfo[1];
    ms->magicPointsProcOwner[0] = maxInfo[0];
    ms->globalIndecesMagicPoints[0] = maxInfo[2];

    rbLifeV::extractSubVector(*(basis[0]), ms->localIndecesMagicPoints,
                              ms->magicPointsProcOwner, interpCoef, 1);
    ms->Qj(0,0) = interpCoef(0);

    VECTOREPETRA rm(*ms->vectorMap);
    VECTOREPETRA QQjrm(*ms->vectorMap);

    Epetra_SerialDenseSolver solverQj;
    solverQj.SetMatrix(ms->Qj);

    for (unsigned int iB = 1; iB < N; iB++)
    {
        rm = *(basis[iB]);
        rm.epetraVector().Update(1., basis[iB]->epetraVector(), 0.);

        computeInterpolationVectorOffline(rm, interpCoef, solverQj, ms);

        computeFeInterpolation(i, j, interpCoef, QQjrm);

        rm.epetraVector().Update(-1., QQjrm.epetraVector(), 1.);

        rbLifeV::getInfNorm(rm, maxInfo);

        ms->magicPointsProcOwner[iB] = maxInfo[0];
        ms->localIndecesMagicPoints[iB] = maxInfo[1];
        ms->globalIndecesMagicPoints[iB] = maxInfo[2];

        ms->Qj.Reshape(iB + 1, iB + 1);

        for (int jB = 0; jB < iB + 1; jB++)
        {
            rbLifeV::extractSubVector(*(basis[jB]), ms->localIndecesMagicPoints,
                                      ms->magicPointsProcOwner, interpCoef, iB + 1);

            for(int kB = 0; kB < iB + 1; kB++)
                ms->Qj(kB, jB) = interpCoef(kB);
        }

        solverQj.SetMatrix(ms->Qj);
    }
}

void
MDEIM::
buildReducedMesh(const unsigned int& i, const unsigned int& j)
{
    identifyReducedNodes(i,j);
    identifyReducedElements(i,j);
}

void
MDEIM::
identifyReducedNodes(const unsigned int& i, const unsigned int& j)
{
    SHP(SingleMDEIMStructure) ms = M_structures(i,j);

    ms->myLocalMagicPoints.resize(ms->N);
    ms->numMyLocalMagicPoints = 0;

    for (int iMp = 0; iMp < ms->N; iMp++)
    {
        if (ms->magicPointsProcOwner[iMp] == M_comm->MyPID())
        {
            ms->myLocalMagicPoints[ms->numMyLocalMagicPoints] = ms->localIndecesMagicPoints[iMp];
            ms->numMyLocalMagicPoints++;
        }
    }

    ms->myLocalMagicPoints.resize(ms->numMyLocalMagicPoints);

    ms->myRowMatrixEntriesOfMagicPoints = new int[ms->numMyLocalMagicPoints];
    ms->myColMatrixEntriesOfMagicPoints = new int[ms->numMyLocalMagicPoints];

    ms->rowLocalReducedIndeces = new int[ms->numMyLocalMagicPoints];

    int localCol(-1);

    auto feMap = M_assembler->getFEspace(i)->map();

    for (int iMp = 0; iMp < ms->numMyLocalMagicPoints; iMp++)
    {
        for (int iR = 0; iR < ms->numMyRows; iR++)
        {
            localCol = ms->myLocalMagicPoints[iMp] - ms->partialSumMyEntries[iR];

            if ((localCol >= 0 ) && (ms->partialSumMyEntries[iR+1] > ms->myLocalMagicPoints[iMp]))
            {
                // attention: we consider the fespace of the row i because that's
                // the range map of the matrix
                ms->myRowMatrixEntriesOfMagicPoints[iMp] = feMap.map(LifeV::Unique)->GID(iR);
                ms->myColMatrixEntriesOfMagicPoints[iMp] = ms->columnIndeces[iR][localCol];

                ms->rowLocalReducedIndeces[iMp] = iR;
            }
        }
    }

    int dimensionOfField = M_assembler->getFEspace(i)->dof().numTotalDof();

    ms->globalReducedNodes = new int[2 * ms->numMyLocalMagicPoints];

    int countMp = 0;

    bool foundCol = false;
    bool foundRow = false;

    for (int iMp = 0; iMp < ms->numMyLocalMagicPoints; iMp++)
    {
        foundRow = false;
        foundCol = (ms->myRowMatrixEntriesOfMagicPoints[iMp] == ms->myColMatrixEntriesOfMagicPoints[iMp]);

        for (int jMp = 0; (jMp < countMp ) && ((!foundRow) || (!foundCol)); jMp++)
        {
            foundRow = (foundRow) || ((ms->globalReducedNodes[jMp]) == (ms->myRowMatrixEntriesOfMagicPoints[iMp] % dimensionOfField));
            foundCol = (foundCol) || ((ms->globalReducedNodes[jMp]) == (ms->myColMatrixEntriesOfMagicPoints[iMp] % dimensionOfField));
        }

        if (!foundRow)
        {
            ms->globalReducedNodes[countMp] = ms->myRowMatrixEntriesOfMagicPoints[iMp] % dimensionOfField;
            countMp++;
        }
        if (!foundCol)
        {
            ms->globalReducedNodes[countMp] = ms->myColMatrixEntriesOfMagicPoints[iMp] % dimensionOfField;
            countMp++;
        }
    }

    ms->numGlobalReducedNodes = countMp;

    // I gather all info by all processors to create a global list of nodes
    int maxGlobalReducedNodes(-1);

    M_comm->MaxAll(&ms->numGlobalReducedNodes, &maxGlobalReducedNodes, 1);

    int* globalReducedNodes = new int[maxGlobalReducedNodes];

    for (int iGrn = 0; iGrn < ms->numGlobalReducedNodes; iGrn++)
    {
        globalReducedNodes[iGrn] = ms->globalReducedNodes[iGrn];
    }
    for (int iGrn = ms->numGlobalReducedNodes; iGrn < maxGlobalReducedNodes; iGrn++)
    {
        globalReducedNodes[iGrn] = -1;
    }

    int numAllGlobalReducedNodes = maxGlobalReducedNodes * M_comm->NumProc();

    int* allGlobalReducedNodes = new int[numAllGlobalReducedNodes];

    M_comm->GatherAll(globalReducedNodes, allGlobalReducedNodes, maxGlobalReducedNodes);

    ms->numMyGlobalReducedNodes = 0;

    for (int iGrn = 0; iGrn < numAllGlobalReducedNodes; iGrn++)
    {
        if (feMap.map(LifeV::Repeated)->MyGID(allGlobalReducedNodes[iGrn]))
            ms->numMyGlobalReducedNodes++;
    }

    int localCount = 0;

    ms->myGlobalReducedNodes = new int[ms->numMyGlobalReducedNodes];

    for (int iGrn = 0; iGrn < numAllGlobalReducedNodes; iGrn++)
    {
        if (feMap.map(LifeV::Repeated)->MyGID(allGlobalReducedNodes[iGrn]))
        {
            ms->myGlobalReducedNodes[localCount] = allGlobalReducedNodes[iGrn];
            localCount++;
        }
    }

    delete [] globalReducedNodes;
    delete [] allGlobalReducedNodes;
}

void
MDEIM::
identifyReducedElements(const unsigned int& i, const unsigned int& j)
{
    SHP(SingleMDEIMStructure) ms = M_structures(i,j);
    auto ufespace = M_assembler->getFEspace(i);

    LifeV::QuadratureRule interpQuad;
    interpQuad.setDimensionShape(shapeDimension(ufespace->refFE().shape()), ufespace->refFE().shape());
    interpQuad.setPoints(ufespace->refFE().refCoor(), std::vector<LifeV::Real>(ufespace->refFE().nbDof(), 0));

    LifeV::CurrentFE interpCFE(ufespace->refFE(), getGeometricMap(*ufespace->mesh()), interpQuad);

    int totalNumberElements = ufespace->mesh()->numElements();
    int numberLocalDof = ufespace->dof().numLocalDof();
    int dimensionOfField = ufespace->dof().numTotalDof();


    ms->numReducedElements = 0;
    bool keepSearching = false;

    // Do the loop over the cells
    for (int iterElement = 0; iterElement < totalNumberElements; iterElement++)
    {
        keepSearching = false;

        // We update the CurrentFE so that we get the coordinates of the nodes
        interpCFE.update(ufespace->mesh()->element(iterElement), LifeV::UPDATE_QUAD_NODES);

        for (int iDim = 0; iDim < ufespace->fieldDim() && !keepSearching; iDim++)
        {
            for (int iterDof = 0; iterDof < numberLocalDof && !keepSearching; iterDof++)
            {
                LifeV::ID globalDofID(ufespace->dof().localToGlobalMap(iterElement, iterDof));
                globalDofID += iDim * dimensionOfField;

                for (int iGrn = 0; iGrn < ms->numMyGlobalReducedNodes && !keepSearching; iGrn++)
                {
                    if (globalDofID == ms->myGlobalReducedNodes[iGrn])
                    {
                        ms->numReducedElements++;
                        keepSearching = true;
                    }
                }
            }
        }
    }

    int localElement = 0;
    ms->reducedElements = new unsigned int [ms->numReducedElements];

    for (int iterElement = 0; iterElement < totalNumberElements; iterElement++)
    {
        interpCFE.update(ufespace->mesh()->element(iterElement), LifeV::UPDATE_QUAD_NODES);

        keepSearching = false;

        for (int iDim = 0; iDim < ufespace->fieldDim() && !keepSearching; iDim++)
        {
            for (int iterDof = 0; iterDof < numberLocalDof  && !keepSearching; iterDof++)
            {
                LifeV::ID globalDofID(ufespace->dof().localToGlobalMap(iterElement, iterDof));

                globalDofID += iDim * dimensionOfField;

                for (int iGrn = 0; iGrn < ms->numMyGlobalReducedNodes  && !keepSearching; iGrn++)
                {
                    if (globalDofID == ms->myGlobalReducedNodes[iGrn])
                    {
                        ms->reducedElements[localElement] = iterElement;
                        localElement++;
                        keepSearching = true;
                    }
                }
            }
        }
    }
}

void
MDEIM::
computeFeInterpolation(const unsigned int& i, const unsigned int& j,
                       Epetra_SerialDenseVector& interpolationCoefficients,
                       VECTOREPETRA& vector)
{
    int currentDeimBasisSize(M_structures(i,j)->Qj.M());
    vector.epetraVector().PutScalar(0.);

    for (int iB = 0; iB < currentDeimBasisSize; iB++)
        vector.epetraVector().Update(interpolationCoefficients(iB),
                                     M_bases(i,j)[iB]->epetraVector(), 1.);
}

void
MDEIM::
computeInterpolationVectorOffline(VECTOREPETRA& vector,
                                  Epetra_SerialDenseVector& interpolationCoefficients,
                                  Epetra_SerialDenseSolver& solver,
                                  SHP(SingleMDEIMStructure) mstruct)
{
    int currentDeimBasisSize(mstruct->Qj.M());

    // Create vector with number of entries to match
    interpolationCoefficients.Resize(currentDeimBasisSize);

    // Extract from real vector the entries to match
    Epetra_SerialDenseVector subVector(currentDeimBasisSize);
    rbLifeV::extractSubVector(vector, mstruct->localIndecesMagicPoints,
                              mstruct->magicPointsProcOwner, subVector, currentDeimBasisSize);

    // Set the unknown vectors and rhs for the interpolation problem
    solver.SetVectors(interpolationCoefficients, subVector);

    // Solve the interpolation coefficients
    solver.Solve( );
}

SHP(VECTOREPETRA)
MDEIM::
vectorizeMatrix(const unsigned int& i, const unsigned int& j, SHP(MATRIXEPETRA) matrix)
{
    SHP(SingleMDEIMStructure) ms = M_structures(i,j);

    SHP(VECTOREPETRA) vectorizedMatrix(new VECTOREPETRA(*ms->vectorMap,
                                                        LifeV::Unique));

    double * values(nullptr);
    int rowCounter = 0;
    int numEntries;

    for (unsigned int iR = 0; iR < ms->numMyRows; iR++)
    {
        values = &(vectorizedMatrix->epetraVector()[0][rowCounter]);
        matrix->matrixPtr()->ExtractMyRowCopy(iR, ms->numMyEntries[iR], numEntries, values);
        rowCounter += ms->numMyEntries[iR];
    }

    vectorizedMatrix->globalAssemble();

    return vectorizedMatrix;
}

void
MDEIM::
vectorizeSnapshots()
{
    M_snapshotsVectorized.resize(M_nRows,M_nCols);
    for (auto snap : M_snapshots)
    {
        for (unsigned int i = 0; i < M_nRows; i++)
        {
            for (unsigned int j = 0; j < M_nCols; j++)
            {
                if (M_structures(i,j))
                    M_snapshotsVectorized(i,j).push_back(vectorizeMatrix(i,j,snap.block(i,j).data()));
            }
        }
    }
}

void
MDEIM::
performPOD()
{
    double podtol = M_data("mdeim/podtol", 1e-5);

    printlog(MAGENTA, "[MDEIM] performing POD(s) ... \n", M_data.getVerbose());

    M_bases.resize(M_nRows, M_nCols);

    for (unsigned int i = 0; i < M_nRows; i++)
    {
        for (unsigned int j = 0; j < M_nCols; j++)
        {
            if (M_structures(i,j))
            {
                auto map = *M_structures(i,j)->vectorMap;
                rbLifeV::ProperOrthogonalDecomposition pod(M_comm,
                                                           map,
                                                           true);
                pod.initPOD(M_snapshotsVectorized(i,j).size(),
                            M_snapshotsVectorized(i,j).data());
                pod.generatePODbasisTol(podtol);

                unsigned int nbfs = pod.getRBdimension();
                std::vector<SHP(VECTOREPETRA)> basisFunctions(nbfs);

                pod.swapReducedBasis(basisFunctions, 0);
                M_bases(i,j) = basisFunctions;

                pickMagicPoints(i,j);
                buildReducedMesh(i,j);
            }
        }
    }

    printlog(MAGENTA, "done\n", M_data.getVerbose());
}

}  // namespace RedMA
