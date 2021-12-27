#include "SparseMatrix.hpp"

namespace RedMA
{

SparseMatrix::
SparseMatrix()
{
}

SparseMatrix::
SparseMatrix(std::vector<shp<DistributedVector>> columnVectors)
{
    const double dropTolerance(2.0 * std::numeric_limits<double>::min());

    unsigned int N = columnVectors.size();
    shp<LifeV::MapEpetra> rangeMap = columnVectors[0]->getVector()->mapPtr();
    shp<Epetra_Comm> comm = rangeMap->commPtr();

    unsigned int myel = N / comm->NumProc();

    // the first process takes care of the remainder
    if (comm->MyPID() == 0)
    {
        myel += N % comm->NumProc();
    }

    shp<LifeV::MapEpetra> domainMap;
    domainMap.reset(new LifeV::MapEpetra(N, myel, 0, comm));

    shp<MATRIXEPETRA> matrix;
    matrix.reset(new MATRIXEPETRA(*rangeMap, N, false));

    Epetra_Map epetraMap = columnVectors[0]->getVector()->epetraMap();
    unsigned int numElements = epetraMap.NumMyElements();

    for (unsigned int i = 0; i < N; i++)
    {
        shp<VECTOREPETRA> clmPtr = columnVectors[i]->getVector();
        VECTOREPETRA columnVectorUnique(*clmPtr, LifeV::Unique);
        for (unsigned int dof = 0; dof < numElements; dof++)
        {
            unsigned int gdof = epetraMap.GID(dof);
            if (columnVectorUnique.isGlobalIDPresent(gdof))
            {
                double value(columnVectorUnique[gdof]);
                if (std::abs(value) > dropTolerance)
                    matrix->addToCoefficient(gdof, i, value);
            }
        }
    }

    comm->Barrier();

    matrix->globalAssemble(domainMap, rangeMap);
    setMatrix(matrix);
}

SparseMatrix::
SparseMatrix(const std::vector<shp<VECTOREPETRA>>& columnVectors)
{
    const double dropTolerance(2.0 * std::numeric_limits<double>::min());

    unsigned int N = columnVectors.size();
    shp<LifeV::MapEpetra> rangeMap = columnVectors[0]->mapPtr();
    shp<Epetra_Comm> comm = rangeMap->commPtr();

    unsigned int myel = N / comm->NumProc();

    // the first process takes care of the remainder
    if (comm->MyPID() == 0)
    {
        myel += N % comm->NumProc();
    }

    shp<LifeV::MapEpetra> domainMap;
    domainMap.reset(new LifeV::MapEpetra(N, myel, 0, comm));

    shp<MATRIXEPETRA> matrix;
    matrix.reset(new MATRIXEPETRA(*rangeMap, N, false));

    Epetra_Map epetraMap = columnVectors[0]->epetraMap();
    unsigned int numElements = epetraMap.NumMyElements();

    for (unsigned int i = 0; i < N; i++)
    {
        VECTOREPETRA columnVectorUnique(*columnVectors[i], LifeV::Unique);
        for (unsigned int dof = 0; dof < numElements; dof++)
        {
            unsigned int gdof = epetraMap.GID(dof);
            if (columnVectorUnique.isGlobalIDPresent(gdof))
            {
                double value(columnVectorUnique[gdof]);
                if (std::abs(value) > dropTolerance)
                {
                    matrix->addToCoefficient(gdof, i, value);
                }
            }
        }
    }

    comm->Barrier();

    matrix->globalAssemble(domainMap, rangeMap);
    setMatrix(matrix);
}

void
SparseMatrix::
add(shp<aMatrix> other)
{
    if (isZero())
    {
        deepCopy(other);
        return;
    }

    if (!other->isZero())
    {
        if (other->nRows() != nRows() || other->nCols() != nCols()) {
            std::string msg = std::string("[SparseMatrix::add] inconsistent dimensions of matrices!\n"
                              "First matrix is ") + std::to_string(nRows()) + std::string("x") +
                              std::to_string(nCols()) + std::string(". Second matrix is ") +
                              std::to_string(other->nRows()) + std::string("x") +
                              std::to_string(other->nCols());
            throw new Exception(msg);
        }

        SparseMatrix* otherMatrix = dynamic_cast<SparseMatrix*>(other.get());
        (*M_matrix) += *otherMatrix->M_matrix;
    }
}

void
SparseMatrix::
multiplyByScalar(const double& coeff)
{
    if (!isZero())
        (*M_matrix) *= coeff;
}

shp<aMatrix>
SparseMatrix::
transpose() const
{
    shp<SparseMatrix> retMatrix(new SparseMatrix());

    if (M_matrix)
    {
        auto domainMap = M_matrix->domainMapPtr();
        auto rangeMap = M_matrix->rangeMapPtr();

        auto transposedMat = M_matrix->transpose();
        transposedMat->globalAssemble(rangeMap, domainMap);

        retMatrix->setMatrix(transposedMat);
    }

    return retMatrix;
}

shp<aMatrix>
SparseMatrix::
multiplyByMatrix(shp<aMatrix> other)
{
    shp<SparseMatrix> retMat(new SparseMatrix());
    if (isZero() || other->isZero())
        return retMat;

    shp<MATRIXEPETRA> otherMatrix = convert<SparseMatrix>(other)->M_matrix;

    if (!isZero() && !other->isZero())
    {
        shp<MATRIXEPETRA> mat(new MATRIXEPETRA(*M_matrix->rangeMapPtr()));
        M_matrix->multiply(false, *otherMatrix, false, *mat, false);
        mat->globalAssemble(otherMatrix->domainMapPtr(), M_matrix->rangeMapPtr());
        retMat->setMatrix(mat);
    }

    return retMat;
}

shp<aVector>
SparseMatrix::
multiplyByVector(shp<aVector> vector)
{
    shp<DistributedVector> vec(new DistributedVector());

    if (isZero())
        return vec;

    if (!isZero() && !vector->isZero())
    {
        if (vector->nRows() != nCols())
            throw new Exception("[SparseMatrix::multiplyByVector] inconsistent dimensions!");

        shp<VECTOREPETRA> otherVector = convert<DistributedVector>(vector)->getVector();
        shp<VECTOREPETRA> res;
        res.reset(new VECTOREPETRA((*M_matrix) * (*otherVector)));
        vec->setVector(res);
    }

    return vec;
}

void
SparseMatrix::
setMatrix(shp<MATRIXEPETRA> matrix)
{
    if (matrix)
    {
        M_matrix = matrix;
        this->M_nRows = M_matrix->rangeMapPtr()->mapSize();
        this->M_nCols = M_matrix->domainMapPtr()->mapSize();
    }
}

bool
SparseMatrix::
isZero() const
{
    return M_matrix == nullptr;
}

void
SparseMatrix::
shallowCopy(shp<aDataWrapper> other)
{
    if (other)
    {
        auto otherMatrix = convert<SparseMatrix>(other);
        setMatrix(otherMatrix->M_matrix);
    }
}

void
SparseMatrix::
deepCopy(shp<aDataWrapper> other)
{
    if (other)
    {
        auto otherMatrix = convert<SparseMatrix>(other);
        shp<MATRIXEPETRA> newMatrix
            (new MATRIXEPETRA(*otherMatrix->M_matrix));
        setMatrix(newMatrix);
    }
}

void
SparseMatrix::
dump(std::string filename) const
{
    M_matrix->spy(filename);
}

shp<DenseMatrix>
SparseMatrix::
toDenseMatrixPtr() const
{
    unsigned int M = M_matrix->rangeMapPtr()->mapSize();
    unsigned int N = M_matrix->domainMapPtr()->mapSize();

    shp<DENSEMATRIX> innerMatrix(new DENSEMATRIX(M,N));

    int maxNumEntries = 0;
    for (unsigned int i = 0; i < M; i++)
        maxNumEntries = std::max(maxNumEntries, M_matrix->matrixPtr()->NumMyEntries(i));

    double* values = new double[maxNumEntries];
    int* indices = new int[maxNumEntries];

    for (unsigned int i = 0; i < M; i++)
    {
        int numMyEntries = M_matrix->matrixPtr()->NumMyEntries(i);
        int numEntries;
        M_matrix->matrixPtr()->ExtractMyRowCopy(i, numMyEntries,
                                                numEntries, values, indices);

        for (unsigned int j = 0; j < numEntries; j++)
            (*innerMatrix)(i,indices[j]) = values[j];
    }

    delete[] values;
    delete[] indices;

    shp<DenseMatrix> retMat(new DenseMatrix());
    retMat->setMatrix(innerMatrix);

    return retMat;
}

DenseMatrix
SparseMatrix::
toDenseMatrix() const
{
    return *toDenseMatrixPtr();
}


SparseMatrix*
SparseMatrix::
clone() const
{
    SparseMatrix* retMatrix = new SparseMatrix();
    if (M_matrix)
    {
        shp<MATRIXEPETRA> newMatrix
            (new MATRIXEPETRA(*M_matrix));
        retMatrix->setMatrix(newMatrix);
    }
    return retMatrix;
}

void
SparseMatrix::
setData(shp<void> data)
{
    setMatrix(spcast<MATRIXEPETRA>(data));
}

shp<void>
SparseMatrix::
data() const
{
    return M_matrix;
}

shp<SparseMatrix>
SparseMatrix::
convertDenseMatrix(shp<DenseMatrix> denseMatrix,
                   shp<Epetra_Comm> comm)
{
    using namespace LifeV;
    shp<SparseMatrix> retMat(new SparseMatrix());

    if (comm->MyPID() != 0)
        throw new Exception("convertDenseVector does not support more than one proc");
    unsigned int rows = denseMatrix->nRows();
    unsigned int cols = denseMatrix->nCols();

    shp<MapEpetra> rangeMap(new MapEpetra(rows, rows, 0, comm));
    shp<MapEpetra> domainMap(new MapEpetra(cols, cols, 0, comm));

    shp<MATRIXEPETRA> matrix(new MATRIXEPETRA(*rangeMap,rows,false));

    for (unsigned int i = 0; i < rows; i++)
    {
        for (unsigned int j = 0; j < cols; j++)
        {
            double value = spcast<DENSEMATRIX>(denseMatrix->data())->operator()(i,j);
            matrix->addToCoefficient(i, j, value);
        }
    }
    matrix->globalAssemble(domainMap, rangeMap);
    retMat->setMatrix(matrix);
    return retMat;
}

shp<MATRIXEPETRA>
SparseMatrix::
getMatrix()
{
    return M_matrix;
}

}
