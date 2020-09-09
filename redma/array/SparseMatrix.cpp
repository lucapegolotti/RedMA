#include "SparseMatrix.hpp"

namespace RedMA
{

SparseMatrix::
SparseMatrix() :
  aMatrix(SPARSE)
{
}

SparseMatrix::
SparseMatrix(const SparseMatrix& other) :
  aMatrix(SPARSE)
{
    if (other.M_matrix)
    {
        // std::shared_ptr<MATRIXEPETRA> newMatrix
        //     (new MATRIXEPETRA(*other.data()));
        setMatrix(other.M_matrix);
    }
}

SparseMatrix::
SparseMatrix(const std::vector<DistributedVector>& columnVectors) :
  aMatrix(SPARSE)
{
    const double dropTolerance(2.0 * std::numeric_limits<double>::min());

    unsigned int N = columnVectors.size();
    std::shared_ptr<LifeV::MapEpetra> rangeMap = static_cast<VECTOREPETRA*>(columnVectors[0].data().get())->mapPtr();
    std::shared_ptr<Epetra_Comm> comm = rangeMap->commPtr();

    unsigned int myel = N / comm->NumProc();

    // the first process takes care of the remainder
    if (comm->MyPID() == 0)
    {
        myel += N % comm->NumProc();
    }

    std::shared_ptr<LifeV::MapEpetra> domainMap;
    domainMap.reset(new LifeV::MapEpetra(N, myel, 0, comm));

    std::shared_ptr<MATRIXEPETRA> matrix;
    matrix.reset(new MATRIXEPETRA(*rangeMap, N, false));

    Epetra_Map epetraMap = static_cast<VECTOREPETRA*>(columnVectors[0].data().get())->epetraMap();
    unsigned int numElements = epetraMap.NumMyElements();

    for (unsigned int i = 0; i < N; i++)
    {
        std::shared_ptr<VECTOREPETRA> clmPtr(static_cast<VECTOREPETRA*>(columnVectors[i].data().get()));
        VECTOREPETRA columnVectorUnique(*clmPtr, LifeV::Unique);
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

SparseMatrix::
SparseMatrix(const std::vector<std::shared_ptr<VECTOREPETRA>>& columnVectors) :
  aMatrix(SPARSE)
{
    const double dropTolerance(2.0 * std::numeric_limits<double>::min());

    unsigned int N = columnVectors.size();
    std::shared_ptr<LifeV::MapEpetra> rangeMap = columnVectors[0]->mapPtr();
    std::shared_ptr<Epetra_Comm> comm = rangeMap->commPtr();

    unsigned int myel = N / comm->NumProc();

    // the first process takes care of the remainder
    if (comm->MyPID() == 0)
    {
        myel += N % comm->NumProc();
    }

    std::shared_ptr<LifeV::MapEpetra> domainMap;
    domainMap.reset(new LifeV::MapEpetra(N, myel, 0, comm));

    std::shared_ptr<MATRIXEPETRA> matrix;
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
add(std::shared_ptr<aMatrix> other)
{
    checkType(other, SPARSE);

    if (isZero())
    {
        hardCopy(other);
        return;
    }

    if (!other->isZero())
    {
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

std::shared_ptr<aMatrix>
SparseMatrix::
transpose() const
{
    std::shared_ptr<SparseMatrix> retMatrix(new SparseMatrix());

    if (!isZero())
    {
        auto domainMap = M_matrix->domainMapPtr();
        auto rangeMap = M_matrix->rangeMapPtr();

        auto transposedMat = M_matrix->transpose();
        transposedMat->globalAssemble(rangeMap, domainMap);

        retMatrix->setMatrix(transposedMat);
    }

    return retMatrix;
}

std::shared_ptr<aMatrix>
SparseMatrix::
multiplyByMatrix(std::shared_ptr<aMatrix> other)
{
    checkType(other, SPARSE);

    std::shared_ptr<MATRIXEPETRA> otherMatrix
    (static_cast<MATRIXEPETRA*>(dynamic_cast<SparseMatrix*>(other.get())->data().get()));

    std::shared_ptr<SparseMatrix> retMat(new SparseMatrix());
    if (!isZero() && other->isZero())
    {
        std::shared_ptr<MATRIXEPETRA> mat(new MATRIXEPETRA(*M_matrix->rangeMapPtr()));
        M_matrix->multiply(false, *otherMatrix, false, *mat, false);
        mat->globalAssemble(otherMatrix->domainMapPtr(), M_matrix->rangeMapPtr());
        retMat->setMatrix(mat);
    }

    return retMat;
}

std::shared_ptr<aVector>
SparseMatrix::
multiplyByVector(std::shared_ptr<aVector> vector)
{
    checkType(vector, DISTRIBUTED);

    std::shared_ptr<DistributedVector> vec(new DistributedVector());
    std::shared_ptr<VECTOREPETRA> otherVector(
    static_cast<VECTOREPETRA*>(dynamic_cast<DistributedVector*>(vector.get())->data().get()));
    if (isZero() && !vector->isZero())
    {
        std::shared_ptr<VECTOREPETRA> res;
        res.reset(new VECTOREPETRA((*M_matrix) * (*otherVector)));

        vec->setVector(res);
    }
    return vec;
}

void
SparseMatrix::
setMatrix(std::shared_ptr<MATRIXEPETRA> matrix)
{
    if (matrix)
    {
        M_matrix = matrix;
        this->M_nRows = M_matrix->rangeMapPtr()->mapSize();
        this->M_nCols = M_matrix->domainMapPtr()->mapSize();
        M_normInf = M_matrix->matrixPtr()->NormInf();
    }
}

bool
SparseMatrix::
isZero() const
{
    if (!M_matrix)
        return true;

    return normInf() < ZEROTHRESHOLD;
}

void
SparseMatrix::
softCopy(std::shared_ptr<aMatrix> other)
{
    if (other)
    {
        checkType(other, SPARSE);
        auto otherMatrix = dynamic_cast<SparseMatrix*>(other.get());
        setMatrix(otherMatrix->M_matrix);
    }
}

void
SparseMatrix::
hardCopy(std::shared_ptr<aMatrix> other)
{
    if (other)
    {
        checkType(other, SPARSE);
        auto otherMatrix = dynamic_cast<SparseMatrix*>(other.get());
        std::shared_ptr<MATRIXEPETRA> newMatrix
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

DenseMatrix
SparseMatrix::
toDenseMatrix() const
{
    unsigned int M = M_matrix->rangeMapPtr()->mapSize();
    unsigned int N = M_matrix->domainMapPtr()->mapSize();

    std::shared_ptr<DENSEMATRIX> innerMatrix(new DENSEMATRIX(M,N));

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

    DenseMatrix retMat;
    retMat.data() = innerMatrix;

    return retMat;
}


SparseMatrix*
SparseMatrix::
clone() const
{
    return new SparseMatrix(*this);
}

std::shared_ptr<void>
SparseMatrix::
data() const
{
    return M_matrix;
}

}
