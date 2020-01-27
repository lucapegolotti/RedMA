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
    std::shared_ptr<LifeV::MapEpetra> domainMap = columnVectors[0].data()->mapPtr();
    std::shared_ptr<Epetra_Comm> comm = domainMap->commPtr();

    unsigned int myel = N / comm->NumProc();

    // the first process takes care of the remainder
    if (comm->MyPID() == 0)
    {
        myel += N % comm->NumProc();
    }

    std::shared_ptr<LifeV::MapEpetra> rangeMap;
    rangeMap.reset(new LifeV::MapEpetra(N, myel, 0, comm));

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

    (*M_matrix) += (*other.data());
    return *this;
}

MatrixEp&
MatrixEp::
operator*=(const double& coeff)
{
    if (!M_matrix)
        (*M_matrix) *= coeff;

    return *this;
}

void
MatrixEp::
hardCopy(const MatrixEp& other)
{
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
    if (!M_matrix)
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

}
