#include "MatrixEp.hpp"

namespace RedMA
{

MatrixEp::
MatrixEp() :
  M_matrix(nullptr)
{
}

MatrixEp::
MatrixEp(const std::vector<VectorEp>& columnVectors)
{
    const double dropTolerance(2.0 * std::numeric_limits<double>::min());

    unsigned int N = columnVectors.size();
    std::shared_ptr<LifeV::MapEpetra> rangeMap = columnVectors[0].data()->mapPtr();
    std::shared_ptr<Epetra_Comm> comm = rangeMap->commPtr();

    unsigned int myel = N / comm->NumProc();

    // the first process takes care of the remainder
    if (comm->MyPID() == 0)
    {
        myel += N % comm->NumProc();
    }

    std::shared_ptr<LifeV::MapEpetra> domainMap;
    domainMap.reset(new LifeV::MapEpetra(N, myel, 0, comm));

    M_matrix.reset(new MATRIXEPETRA(*rangeMap, N, false));

    Epetra_Map epetraMap = columnVectors[0].data()->epetraMap();
    unsigned int numElements = epetraMap.NumMyElements();

    for (unsigned int i = 0; i < N; i++)
    {
        VECTOREPETRA columnVectorUnique(*columnVectors[i].data(), LifeV::Unique);
        for (unsigned int dof = 0; dof < numElements; dof++)
        {
            unsigned int gdof = epetraMap.GID(dof);
            if (columnVectorUnique.isGlobalIDPresent(gdof))
            {
                double value(columnVectorUnique[gdof]);
                if (std::abs(value) > dropTolerance)
                {
                    M_matrix->addToCoefficient(gdof, i, value);
                }
            }
        }
    }

    comm->Barrier();

    M_matrix->globalAssemble(domainMap, rangeMap);
}

MatrixEp::
MatrixEp(const std::vector<std::shared_ptr<VECTOREPETRA>>& columnVectors)
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

    M_matrix.reset(new MATRIXEPETRA(*rangeMap, N, false));

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
                    M_matrix->addToCoefficient(gdof, i, value);
                }
            }
        }
    }

    comm->Barrier();

    M_matrix->globalAssemble(domainMap, rangeMap);
}

MatrixEp
MatrixEp::
transpose() const
{
    MatrixEp retMatrix;

    if (M_matrix)
    {
        auto domainMap = M_matrix->domainMapPtr();
        auto rangeMap = M_matrix->rangeMapPtr();

        retMatrix.data() = M_matrix->transpose();
        retMatrix.data()->globalAssemble(rangeMap, domainMap);
    }

    return retMatrix;
}

MatrixEp
MatrixEp::
operator*(const MatrixEp& other)
{
    MatrixEp mat;

    if (M_matrix && other.data())
    {
        mat.data().reset(new MATRIXEPETRA(*M_matrix->rangeMapPtr()));
        M_matrix->multiply(false, *other.data(), false, *mat.data(), false);
        mat.data()->globalAssemble(other.data()->domainMapPtr(), M_matrix->rangeMapPtr());
    }

    return mat;
}

MatrixEp
MatrixEp::
operator+(const MatrixEp& other)
{
    MatrixEp retMatrix;

    if (!M_matrix)
    {
        hardCopy(other);
        return other;
    }

    retMatrix.data().reset(new MATRIXEPETRA(*M_matrix));

    (*retMatrix.data()) += (*other.data());
    return retMatrix;
}

MatrixEp&
MatrixEp::
operator+=(const MatrixEp& other)
{
    if (!M_matrix)
    {
        hardCopy(other);
        return *this;
    }

    if (other.data())
        (*M_matrix) += (*other.data());
    return *this;
}

MatrixEp&
MatrixEp::
operator-=(const MatrixEp& other)
{
    if (!M_matrix)
    {
        hardCopy(other);
        *this *= (-1.0);
        return *this;
    }

    if (other.data())
        (*M_matrix) -= (*other.data());
    return *this;
}

MatrixEp&
MatrixEp::
operator*=(const double& coeff)
{
    if (M_matrix)
        (*M_matrix) *= coeff;

    return *this;
}

void
MatrixEp::
hardCopy(const MatrixEp& other)
{
    if (other.data())
        M_matrix.reset(new MATRIXEPETRA(*other.data()));
}

void
MatrixEp::
softCopy(const MatrixEp& other)
{
    M_matrix = other.M_matrix;
}

VectorEp
MatrixEp::
operator*(const VectorEp& vector)
{
    VectorEp vec;
    if (M_matrix && vector.data())
    {
        std::shared_ptr<VECTOREPETRA> res;
        res.reset(new VECTOREPETRA((*M_matrix) * (*vector.data())));

        vec.data() = res;
    }
    return vec;
}

void
MatrixEp::
getRowProperty(std::shared_ptr<LifeV::MapEpetra>& outMap)
{
    outMap.reset(new LifeV::MapEpetra(M_matrix->domainMap()));
}

void
MatrixEp::
getColProperty(std::shared_ptr<LifeV::MapEpetra>& outMap)
{
    outMap.reset(new LifeV::MapEpetra(M_matrix->rangeMap()));
}

std::shared_ptr<MATRIXEPETRA>&
MatrixEp::
data()
{
    return M_matrix;
}

std::shared_ptr<MATRIXEPETRA>
MatrixEp::
data() const
{
    return M_matrix;
}

void
MatrixEp::
dump(std::string filename) const
{
    M_matrix->spy(filename);
}

DenseMatrix
MatrixEp::
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

}
