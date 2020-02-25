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
            }
        }
    }

    printlog(MAGENTA, "done\n", M_data.getVerbose());
}

}  // namespace RedMA
