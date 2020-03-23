#include "MDEIM.hpp"

namespace RedMA
{

MDEIM::
MDEIM() :
  M_isInitialized(false)
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
addSnapshot(MatrixEp newSnapshot)
{
    M_snapshots.push_back(newSnapshot);
}

void
MDEIM::
setFESpace(SHP(FESPACE) fespace)
{
    M_fespace = fespace;
}

void
MDEIM::
initialize(MatrixEp matrix)
{
    SHP(MATRIXEPETRA) matrixEpetra = matrix.data();

    if (matrixEpetra)
    {
        SHP(MDEIMStructure) ms(new MDEIMStructure());

        ms->numGlobalNonzeros = matrixEpetra->matrixPtr()->NumGlobalNonzeros();
        ms->numMyNonzeros = matrixEpetra->matrixPtr()->NumMyNonzeros();
        ms->numMyRows = matrixEpetra->matrixPtr()->NumMyRows();

        ms->numMyEntries.resize(ms->numMyRows);
        ms->columnIndices.resize(ms->numMyRows);
        int numEntries;

        for (unsigned int iR = 0; iR < ms->numMyRows; iR++)
        {
            ms->numMyEntries[iR] = matrixEpetra->matrixPtr()->NumMyEntries(iR);
            ms->columnIndices[iR].resize(ms->numMyEntries[iR]);
            double * values = new double[ms->numMyEntries[iR]];
            matrixEpetra->matrixPtr()->ExtractMyRowCopy(iR, ms->numMyEntries[iR],
                                                        numEntries, values,
                                                        ms->columnIndices[iR].data());

            for(unsigned int iC; iC < numEntries; iC++)
            {
                ms->columnIndices[iR][iC] = matrixEpetra->matrixPtr()->ColMap().
                                                    GID(ms->columnIndices[iR][iC]);
            }
            delete [] values;
        }

        ms->partialSumMyEntries.resize(ms->numMyRows + 1);
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
        M_isInitialized = true;
        M_structure = ms;
    }
}
//
void
MDEIM::
performMDEIM(std::string outdir)
{
    if (M_isInitialized)
    {
        vectorizeSnapshots();
        performPOD(outdir);
    }
}

void
MDEIM::
pickMagicPoints()
{
    SHP(MDEIMStructure) ms = M_structure;

    // get first coefficient manually
    int maxInfo[3];
    Epetra_SerialDenseVector interpCoef(1);

    unsigned int N = M_basis.size();

    ms->N = N;
    ms->Qj.Reshape(1,1);

    // I keep these structures even if I don't plan to call this routine in parallel
    ms->localIndicesMagicPoints.resize(N);
    ms->magicPointsProcOwner.resize(N);
    ms->globalIndicesMagicPoints.resize(N);

    rbLifeV::getInfNorm(*(M_basis[0]), maxInfo);

    ms->localIndicesMagicPoints[0] = maxInfo[1];
    ms->magicPointsProcOwner[0] = maxInfo[0];
    ms->globalIndicesMagicPoints[0] = maxInfo[2];

    rbLifeV::extractSubVector(*(M_basis[0]), ms->localIndicesMagicPoints,
                              ms->magicPointsProcOwner, interpCoef, 1);
    ms->Qj(0,0) = interpCoef(0);

    VECTOREPETRA rm(*ms->vectorMap);
    VECTOREPETRA QQjrm(*ms->vectorMap);

    Epetra_SerialDenseSolver solverQj;
    solverQj.SetMatrix(ms->Qj);

    for (unsigned int iB = 1; iB < N; iB++)
    {
        rm = *(M_basis[iB]);
        rm.epetraVector().Update(1., M_basis[iB]->epetraVector(), 0.);

        computeInterpolationVectorOffline(rm, interpCoef, solverQj);
        std::cout << "interpCoef = " << interpCoef.Norm2() << std::endl << std::flush;
        computeFeInterpolation(interpCoef, QQjrm);
        std::cout << "QQjrm = " << QQjrm.norm2() << std::endl << std::flush;

        rm.epetraVector().Update(-1., QQjrm.epetraVector(), 1.);

        rbLifeV::getInfNorm(rm, maxInfo);

        ms->magicPointsProcOwner[iB] = maxInfo[0];
        ms->localIndicesMagicPoints[iB] = maxInfo[1];
        ms->globalIndicesMagicPoints[iB] = maxInfo[2];

        ms->Qj.Reshape(iB + 1, iB + 1);

        for (int jB = 0; jB < iB + 1; jB++)
        {
            // std::cout << jB << "=========" << std::endl << std::flush;
            // rbLifeV::extractSubVector(*(M_basis[jB]), ms->localIndicesMagicPoints,
            //                           ms->magicPointsProcOwner, interpCoef, iB + 1);

            interpCoef.Resize(iB + 1);

            for (int iV = 0; iV < iB + 1; iV++)
            {
                interpCoef[iV] = M_basis[jB]->epetraVector().Values()[ms->localIndicesMagicPoints[iV]];
                M_basis[jB]->mapPtr()->commPtr()->Broadcast(&(interpCoef[iV]), 1, ms->magicPointsProcOwner[iV]);

                if (std::abs(interpCoef[iV]) < 1e-15)
                {
                    std::cout << "basis value " << M_basis[jB]->epetraVector().Values()[ms->localIndicesMagicPoints[iV]] << std::endl << std::flush;
                    std::cout << "approximation value " << QQjrm.epetraVector().Values()[ms->localIndicesMagicPoints[iV]] << std::endl << std::flush;
                }
            }

            // for (int kB = 0; kB < iB + 1; kB++)
            //     std::cout << "vector = " << interpCoef(kB) << std::endl << std::flush;

            for(int kB = 0; kB < iB + 1; kB++)
                ms->Qj(kB, jB) = interpCoef(kB);
        }
        std::cout << "Qj norm = " << ms->Qj.NormInf() << std::endl << std::flush;

        solverQj.SetMatrix(ms->Qj);
    }

    unsigned int M = ms->Qj.M();
    N = ms->Qj.N();

    std::ostringstream streamObj;
    for (unsigned int i = 0; i < M; i++)
    {
        for (unsigned int j = 0; j < N; j++)
        {
            streamObj << M_structure->Qj(i,j);
            if (j != M-1)
                streamObj << ",";
        }
        if (i != N-1)
            streamObj << "\n";
    }
    std::cout << streamObj.str() << std::endl << std::flush;
}

void
MDEIM::
buildReducedMesh()
{
    identifyReducedNodes();
    identifyReducedElements();
}

void
MDEIM::
prepareOnline()
{
    if (M_isInitialized)
    {
        SHP(MDEIMStructure) ms = M_structure;
        auto matrixEpetra = M_snapshots[0].data();

        buildReducedMesh();

        ms->numReducedGlobalNonzeros = matrixEpetra->matrixPtr()->NumGlobalNonzeros();
        ms->numReducedMyNonzeros = matrixEpetra->matrixPtr()->NumMyNonzeros();
        ms->numReducedMyRows = matrixEpetra->matrixPtr()->NumMyRows();

        ms->numMyReducedEntries.resize(ms->numMyLocalMagicPoints);
        ms->columnLocalReducedIndices.resize(ms->numMyLocalMagicPoints);
        int numEntries;

        for (int iMp = 0; iMp < ms->numMyLocalMagicPoints; iMp++)
        {
            ms->numMyReducedEntries[iMp] = matrixEpetra->matrixPtr()->NumMyEntries(ms->rowLocalReducedIndices[iMp]);
            double* values = new double[ms->numMyReducedEntries[iMp]];
            int* indeces = new int[ms->numMyReducedEntries[iMp]];

            matrixEpetra->matrixPtr()->ExtractMyRowCopy(ms->rowLocalReducedIndices[iMp], ms->numMyReducedEntries[iMp], numEntries, values, indeces);

            for (unsigned int iC = 0; iC < ms->numMyReducedEntries[iMp]; iC++)
            {
                if (ms->myColMatrixEntriesOfMagicPoints[iMp] == matrixEpetra->matrixPtr()->ColMap().GID(indeces[iC]))
                {
                    ms->columnLocalReducedIndices[iMp] = iC;
                    iC = ms->numMyReducedEntries[iMp];
                }
            }
            delete [] values;
            delete [] indeces;
        }
    }
}

void
MDEIM::
identifyReducedNodes()
{
    SHP(MDEIMStructure) ms = M_structure;

    ms->myLocalMagicPoints.resize(ms->N);
    ms->numMyLocalMagicPoints = 0;

    for (int iMp = 0; iMp < ms->N; iMp++)
    {
        if (ms->magicPointsProcOwner[iMp] == M_comm->MyPID())
        {
            ms->myLocalMagicPoints[ms->numMyLocalMagicPoints] = ms->localIndicesMagicPoints[iMp];
            ms->numMyLocalMagicPoints++;
        }
    }

    ms->myLocalMagicPoints.resize(ms->numMyLocalMagicPoints);

    ms->myRowMatrixEntriesOfMagicPoints.resize(ms->numMyLocalMagicPoints);
    ms->myColMatrixEntriesOfMagicPoints.resize(ms->numMyLocalMagicPoints);

    ms->rowLocalReducedIndices.resize(ms->numMyLocalMagicPoints);

    int localCol = -1;

    auto feMap = M_fespace->map();

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
                ms->myColMatrixEntriesOfMagicPoints[iMp] = ms->columnIndices[iR][localCol];

                ms->rowLocalReducedIndices[iMp] = iR;
            }
        }
    }

    int dimensionOfField = M_fespace->dof().numTotalDof();

    ms->globalReducedNodes.resize(2 * ms->numMyLocalMagicPoints);

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

    ms->myGlobalReducedNodes.resize(ms->numMyGlobalReducedNodes);

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
identifyReducedElements()
{
    SHP(MDEIMStructure) ms = M_structure;

    auto ufespace = M_fespace;

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
    ms->reducedElements.resize(ms->numReducedElements);

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
computeInterpolationVectorOnline(Epetra_SerialDenseVector& interpVector,
                                 MatrixEp reducedMat)
{
    std::cout << "+++++++++" << std::endl << std::flush;
    auto ms = M_structure;

    Epetra_SerialDenseVector rhsVector;

    computeInterpolationRhsOnline(rhsVector, reducedMat);
    std::cout << "rhsVector " << rhsVector.Norm2() << std::endl << std::flush;
    Epetra_SerialDenseSolver solverQj;

    // for some reason the matrix is changed after solve. So, we make a copy of Qj
    Epetra_SerialDenseMatrix Qj = M_structure->Qj;

    solverQj.SetMatrix(Qj);
    std::cout << "Qj norm = " << Qj.NormInf() << std::endl << std::flush;
    solverQj.SetVectors(interpVector, rhsVector);
    solverQj.Solve();
    std::cout << "interpVector " << interpVector.Norm2() << std::endl << std::flush;

    if (interpVector.Norm2() != interpVector.Norm2())
    {
        std::cout << "it's nan, quitting" << std::endl << std::flush;

        std::ofstream outfile;
        outfile.open("Qj.csv", std::ios_base::out);

        unsigned int M = Qj.M();
        unsigned int N = Qj.N();

        std::ostringstream streamObj;
        for (unsigned int i = 0; i < M; i++)
        {
            for (unsigned int j = 0; j < N; j++)
            {
                streamObj << M_structure->Qj(i,j);
                if (j != M-1)
                    streamObj << ",";
            }
            if (i != N-1)
                streamObj << "\n";
        }
        outfile << streamObj.str();
        outfile.close();

        outfile.open("rhsVector.csv", std::ios_base::out);
        N = rhsVector.Length();

        streamObj.str() = "";
        for (unsigned int i = 0; i < N; i++)
        {
            streamObj << rhsVector(i);
            if (i != N-1)
                streamObj << "\n";
        }
        outfile << streamObj.str();
        outfile.close();

        exit(1);
    }
}

void
MDEIM::
computeInterpolationRhsOnline(Epetra_SerialDenseVector& interpVector,
                              MatrixEp reducedMat)
{
    auto ms = M_structure;

    interpVector.Resize(ms->N);

    auto map = M_fespace->map();

    SHP(MATRIXEPETRA) Ah = reducedMat.data();

    Epetra_SerialDenseVector localInterpVector(std::max(ms->numMyLocalMagicPoints, 1));
    localInterpVector(0) = 1.2345;

    int numEntries;
    double * values;

    for (int iMp = 0; iMp < ms->numMyLocalMagicPoints; iMp++)
    {
        values = new double[ms->numMyReducedEntries[iMp]];

        Ah->matrixPtr()->ExtractMyRowCopy(ms->rowLocalReducedIndices[iMp],
                                          ms->numMyReducedEntries[iMp],
                                          numEntries, values);
        localInterpVector(iMp) = values[ms->columnLocalReducedIndices[iMp]];
        delete[] values;
    }

    unsigned int moveOn = 0;
    double value;

    for (int iMp = 0; iMp < ms->N; iMp++)
    {
        value = localInterpVector(moveOn);
        M_comm->Broadcast(&value, 1, ms->magicPointsProcOwner[iMp]);

        interpVector(iMp) = value;

        if (ms->magicPointsProcOwner[iMp] == M_comm->MyPID())
            moveOn++;
    }
}

void
MDEIM::
computeFeInterpolation(DENSEVECTOR& interpolationCoefficients,
                       VECTOREPETRA& vector)
{
    int currentDeimBasisSize(M_structure->Qj.M());
    vector.epetraVector().PutScalar(0.);

    for (int iB = 0; iB < currentDeimBasisSize; iB++)
        vector.epetraVector().Update(interpolationCoefficients(iB),
                                     M_basis[iB]->epetraVector(), 1.);
}

void
MDEIM::
computeProjectedInterpolation(DENSEVECTOR& interpolationCoefficients,
                              DENSEVECTOR& vector)
{
    int currentDeimBasisSize(M_structure->Qj.M());
    vector.Scale(0.);

    for (int iB = 0; iB < currentDeimBasisSize; iB++)
    {
        DENSEVECTOR aux(*M_basisProjected[iB]);
        aux.Scale(interpolationCoefficients(iB));
        vector += aux;
    }
}

void
MDEIM::
checkOnline(MatrixEp reducedMatrix, MatrixEp fullMatrix)
{
    if (M_isInitialized)
    {
        FEMATRIX approxMat(assembleMatrix(reducedMatrix));

        SHP(MATRIXEPETRA) actualMatrix = fullMatrix.data();

        std::cout << "=======CHECKING ONLINE MDEIM========" << std::endl << std::flush;
        std::cout << "NORM apprMatrix = " << approxMat.data()->normFrobenius() << std::endl << std::flush;

        approxMat -= fullMatrix;

        std::cout << "NORM actualMatrix = " << actualMatrix->normFrobenius() << std::endl << std::flush;
        std::cout << "REL NORM DIFFERENCE = " << approxMat.data()->normFrobenius() / actualMatrix->normFrobenius() << std::endl << std::flush;
    }
}

FEMATRIX
MDEIM::
assembleMatrix(FEMATRIX reducedMatrix)
{
    FEMATRIX retMat;

    if (M_isInitialized)
    {
        auto ms = M_structure;

        // Comparison vectors
        SHP(VECTOREPETRA) approximation;

        // Compute interpolation vector
        Epetra_SerialDenseVector myInterpVector(ms->N);
        computeInterpolationVectorOnline(myInterpVector, reducedMatrix);
        std::cout << "reducedMatrix " << reducedMatrix.data()->normInf() << std::endl << std::flush;
        // Build FEM vector from interpolation vector
        approximation.reset(new VECTOREPETRA(*ms->vectorMap));
        computeFeInterpolation(myInterpVector, *approximation);
        std::cout << "myInterpVector " << myInterpVector.Norm2() << std::endl << std::flush;
        std::cout << "approximation " << approximation->norm2() << std::endl << std::flush;
        SHP(MATRIXEPETRA) apprMatrix;
        apprMatrix.reset(new MATRIXEPETRA(M_fespace->map(), 100));
        apprMatrix->matrixPtr( )->Scale(0.);

        reconstructMatrixFromVectorizedForm(*approximation, *apprMatrix);
        apprMatrix->globalAssemble(M_domainMap, M_rangeMap);
        std::cout << "apprMatrix " << apprMatrix->normInf() << std::endl << std::flush;
        retMat.data() = apprMatrix;
    }
    return retMat;
}

RBMATRIX
MDEIM::
assembleProjectedMatrix(FEMATRIX reducedMatrix)
{
    RBMATRIX retMat;

    if (M_isInitialized)
    {
        auto ms = M_structure;

        // Comparison vectors
        SHP(DENSEVECTOR) approximation;

        // Compute interpolation vector
        Epetra_SerialDenseVector myInterpVector(ms->N);
        computeInterpolationVectorOnline(myInterpVector, reducedMatrix);
        // Build FEM vector from interpolation vector
        approximation.reset(new DENSEVECTOR((ms->Nleft) * (ms->Nright)));
        computeProjectedInterpolation(myInterpVector, *approximation);

        SHP(DENSEMATRIX) apprMatrix(new DENSEMATRIX(ms->Nleft, ms->Nright));
        apprMatrix->Scale(0.);

        reconstructMatrixFromVectorizedForm(*approximation, *apprMatrix);

        retMat.data() = apprMatrix;
    }
    return retMat;
}

void
MDEIM::
loadMDEIM(std::string pathdir)
{
    using namespace boost::filesystem;

    if (exists(pathdir + "/structure.mstr"))
    {
        M_structure.reset(new MDEIMStructure(pathdir + "/structure.mstr", M_comm));

        // create_directory(pathdir + "/copy");
        // M_structure->dumpMDEIMStructure(pathdir + "/copy");
        M_isInitialized = true;
    }

    if (exists(pathdir + "/basis.mbasis") && M_data("rb/online/mdeim/loadfullbasis", false))
    {
        if (!M_isInitialized)
            throw new Exception("MDEIM structure not loaded!");
        loadBasis(pathdir + "/basis.mbasis");
    }

    if (exists(pathdir + "/projbasis.mbasis"))
    {
        if (!M_isInitialized)
            throw new Exception("MDEIM structure not loaded!");
        loadProjectedBasis(pathdir + "/projbasis.mbasis");
    }

}

void
MDEIM::
loadBasis(std::string filename)
{
    std::ifstream infile(filename);
    std::string line;
    while(std::getline(infile,line))
    {
        SHP(VECTOREPETRA) newVector(new VECTOREPETRA(*M_structure->vectorMap));

        std::stringstream linestream(line);
        std::string value;
        unsigned int i = 0;
        while(getline(linestream,value,','))
        {
            newVector->operator[](i) = std::atof(value.c_str());
            i++;
        }
        if (i != newVector->epetraVector().GlobalLength())
            throw new Exception("Stored snapshot length does not match fespace dimension!");
        M_basis.push_back(newVector);
    }
    infile.close();
}

void
MDEIM::
loadProjectedBasis(std::string filename)
{
    std::ifstream infile(filename);
    std::string line;
    while(std::getline(infile,line))
    {
        unsigned int N = (M_structure->Nleft) * (M_structure->Nright);
        SHP(DENSEVECTOR) newVector(new DENSEVECTOR(N));
        std::stringstream linestream(line);
        std::string value;
        unsigned int i = 0;
        while(getline(linestream,value,','))
        {
            newVector->operator()(i) = std::atof(value.c_str());
            i++;
        }
        if (i != N)
            throw new Exception("Stored snapshot length does not match fespace dimension!");
        M_basisProjected.push_back(newVector);
    }
    infile.close();
}

void
MDEIM::
reconstructMatrixFromVectorizedForm(VECTOREPETRA& vectorizedAh,
                                    MATRIXEPETRA& Ah)
{
    auto ms = M_structure;

    int rowStart = 0;

    for (int iR = 0; iR < ms->numMyRows; iR++)
    {
        Ah.matrixPtr()->InsertGlobalValues(Ah.matrixPtr()->GRID(iR), ms->numMyEntries[iR],
                                           vectorizedAh.epetraVector()[0] + rowStart, ms->columnIndices[iR].data());
        rowStart += ms->numMyEntries[iR];
    }
}

void
MDEIM::
reconstructMatrixFromVectorizedForm(DENSEVECTOR& vectorizedAn,
                                    DENSEMATRIX& An)
{
    auto ms = M_structure;

    unsigned int nrows = ms->Nleft;
    unsigned int ncols = ms->Nright;

    An.Reshape(nrows, ncols);

    unsigned int linearcount = 0;
    for (unsigned int i = 0; i < nrows; i++)
    {
        for (unsigned int j = 0; j < ncols; j++)
        {
            An(i,j) = vectorizedAn(linearcount);
            linearcount++;
        }
    }
}

void
MDEIM::
computeInterpolationVectorOffline(VECTOREPETRA& vector,
                                  Epetra_SerialDenseVector& interpolationCoefficients,
                                  Epetra_SerialDenseSolver& solver)
{
    int currentDeimBasisSize(M_structure->Qj.M());

    // Create vector with number of entries to match
    interpolationCoefficients.Resize(currentDeimBasisSize);

    // Extract from real vector the entries to match
    Epetra_SerialDenseVector subVector(currentDeimBasisSize);
    rbLifeV::extractSubVector(vector, M_structure->localIndicesMagicPoints,
                              M_structure->magicPointsProcOwner, subVector,
                              currentDeimBasisSize);

    // Set the unknown vectors and rhs for the interpolation problem
    solver.SetVectors(interpolationCoefficients, subVector);

    auto Qj = solver.Matrix();
    // Solve the interpolation coefficients
    solver.Solve( );

    if (interpolationCoefficients.Norm2() != interpolationCoefficients.Norm2())
    {

        std::cout << "it's nan, quitting" << std::endl << std::flush;

        std::ofstream outfile;
        outfile.open("matrixInterp.csv", std::ios_base::out);

        unsigned int M = Qj->M();
        unsigned int N = Qj->N();

        std::ostringstream streamObj;
        for (unsigned int i = 0; i < M; i++)
        {
            for (unsigned int j = 0; j < N; j++)
            {
                streamObj << Qj->operator()(i,j);
                if (j != M-1)
                    streamObj << ",";
            }
            if (i != N-1)
                streamObj << "\n";
        }
        outfile << streamObj.str();
        outfile.close();

        outfile.open("subVector.csv", std::ios_base::out);
        N = subVector.Length();

        streamObj.str() = "";
        for (unsigned int i = 0; i < N; i++)
        {
            streamObj << subVector(i);
            if (i != N-1)
                streamObj << "\n";
        }
        outfile << streamObj.str();
        outfile.close();

        exit(1);
    }
}
//
SHP(VECTOREPETRA)
MDEIM::
vectorizeMatrix(MatrixEp matrix)
{
    SHP(VECTOREPETRA) vectorizedMatrix(new VECTOREPETRA(*M_structure->vectorMap,
                                                        LifeV::Unique));

    double* values = nullptr;
    int rowCounter = 0;
    int numEntries;

    for (unsigned int iR = 0; iR < M_structure->numMyRows; iR++)
    {
        values = &(vectorizedMatrix->epetraVector()[0][rowCounter]);
        matrix.data()->matrixPtr()->ExtractMyRowCopy(iR, M_structure->numMyEntries[iR],
                                                     numEntries, values);
        rowCounter += M_structure->numMyEntries[iR];
    }

    vectorizedMatrix->globalAssemble();

    return vectorizedMatrix;
}

void
MDEIM::
vectorizeSnapshots()
{
    for (auto snap : M_snapshots)
        M_snapshotsVectorized.push_back(vectorizeMatrix(snap));
}

void
MDEIM::
performPOD(std::string outdir)
{
    double podtol = M_data("rb/offline/mdeim/podtol", 1e-5);

    printlog(MAGENTA, "[MDEIM] performing POD(s) ... \n", M_data.getVerbose());

    if (M_isInitialized)
    {
        auto map = *M_structure->vectorMap;
        rbLifeV::ProperOrthogonalDecomposition pod(M_comm, map, true);
        pod.initPOD(M_snapshotsVectorized.size(),
                    M_snapshotsVectorized.data());
        pod.setSvdFileName(outdir + "/svd.txt");
        pod.generatePODbasisTol(podtol);

        unsigned int nbfs = pod.getRBdimension();
        std::vector<SHP(VECTOREPETRA)> basisFunctions(nbfs);

        pod.swapReducedBasis(basisFunctions, 0);
        M_basis = basisFunctions;

        pickMagicPoints();
        buildReducedMesh();
    }

    printlog(MAGENTA, "done\n", M_data.getVerbose());
}

void
MDEIM::
checkOnSnapshots()
{
    for (auto snap : M_snapshots)
        checkOnline(snap, snap);
}

void
MDEIM::
dumpMDEIM(std::string dir)
{
    if (M_isInitialized)
    {
        M_structure->dumpMDEIMStructure(dir);
        dumpBasis(dir);
        if (M_basisProjected.size() > 0)
            dumpProjectedBasis(dir);
    }
}

void
MDEIM::
dumpBasis(std::string dir)
{
    std::ofstream outfile;
    outfile.open(dir + "/basis.mbasis", std::ios_base::out);

    for (auto basis : M_basis)
    {
        VectorEp newVec;
        newVec.data() = basis;
        outfile << newVec.getString(',') + "\n";
    }

    outfile.close();
}

void
MDEIM::
dumpProjectedBasis(std::string dir)
{
    std::ofstream outfile;
    outfile.open(dir + "/projbasis.mbasis", std::ios_base::out);

    for (auto basis : M_basisProjected)
    {
        DenseVector newVec;
        newVec.data() = basis;
        outfile << newVec.getString(',') + "\n";
    }

    outfile.close();
}

void
MDEIM::
projectMDEIM(std::vector<SHP(VECTOREPETRA)> leftBasis,
             std::vector<SHP(VECTOREPETRA)> rightBasis)
{
    if (M_isInitialized)
    {
        unsigned int Nleft = leftBasis.size();
        unsigned int Nright = rightBasis.size();

        M_structure->Nleft = Nleft;
        M_structure->Nright = Nright;

        for (unsigned int bind = 0; bind < M_basis.size(); bind++)
        {
            SHP(DENSEMATRIX) newMatrix(new DENSEMATRIX(Nleft, Nright));

            SHP(MATRIXEPETRA) basisMatrix(new MATRIXEPETRA(*M_rangeMap));

            reconstructMatrixFromVectorizedForm(*M_basis[bind], *basisMatrix);
            basisMatrix->globalAssemble(M_domainMap, M_rangeMap);
            for (unsigned int i = 0; i < Nleft; i++)
            {
                VECTOREPETRA aux(*rightBasis[0]->mapPtr());
                basisMatrix->matrixPtr()->Multiply(true,
                                leftBasis[i]->epetraVector(), aux.epetraVector());

                for (unsigned int j = 0; j < Nright; j++)
                {
                    (*newMatrix)(i,j) += aux.dot(*rightBasis[j]);
                }
            }

            SHP(DENSEVECTOR) newBasis(new DENSEVECTOR(Nleft * Nright));

            unsigned int count = 0;
            for (unsigned int i = 0; i < Nleft; i++)
            {
                for (unsigned int j = 0; j < Nright; j++)
                {
                    (*newBasis)(count) = (*newMatrix)(i,j);
                    count++;
                }
            }
            M_basisProjected.push_back(newBasis);
        }
    }
}

}  // namespace RedMA
