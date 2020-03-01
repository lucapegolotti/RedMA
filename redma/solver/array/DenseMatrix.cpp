#include "DenseMatrix.hpp"

namespace RedMA
{

DenseMatrix::
DenseMatrix() :
  M_matrix(nullptr)
{
}

// DenseMatrix::
// DenseMatrix(const std::vector<int>& columnVectors)
// {
//     auto vecs = columnVectors;
//     // unsigned int ncols = columnVectors[0].data()->Length();
//     //
//     // M_matrix->Reshape(nrows,ncols);
//     //
//     // for (unsigned int i = 0; i < nrows; i++)
//     // {
//     //     for (unsigned int j = 0; j < ncols; j++)
//     //     {
//     //         (*matrix)(i,j) = columnVectors[i][j];
//     //     }
//     // }
// }

DenseMatrix
DenseMatrix::
transpose() const
{
    DenseMatrix retMatrix;

    if (M_matrix)
    {
        unsigned int nrows = M_matrix->N();
        unsigned int ncols = M_matrix->M();

        retMatrix.data()->Reshape(nrows,ncols);

        for (unsigned int i = 0; i < nrows; i++)
        {
            for (unsigned int j = 0; j < ncols; j++)
            {
                (*retMatrix.M_matrix)(i,j) = (*M_matrix)(j,i);
            }
        }
    }

    return retMatrix;
}

DenseMatrix
DenseMatrix::
operator*(const DenseMatrix& other)
{
    DenseMatrix mat;

    if (M_matrix && other.data())
    {
        mat.data().reset(new DENSEMATRIX(M_matrix->M(), other.M_matrix->N()));
        mat.M_matrix->Multiply('N', 'N', 1.0, *M_matrix, *other.M_matrix, 0.0);
    }

    return mat;
}

DenseMatrix
DenseMatrix::
operator+(const DenseMatrix& other)
{
    DenseMatrix retMatrix;

    if (!M_matrix)
    {
        hardCopy(other);
        return other;
    }

    retMatrix.data().reset(new DENSEMATRIX(*M_matrix));

    (*retMatrix.data()) += (*other.data());
    return retMatrix;
}

DenseMatrix&
DenseMatrix::
operator+=(const DenseMatrix& other)
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

DenseMatrix&
DenseMatrix::
operator-=(const DenseMatrix& other)
{
    if (!M_matrix)
    {
        hardCopy(other);
        *this *= (-1.0);
        return *this;
    }

    if (other.data())
    {
        DENSEMATRIX otherCopy(*other.data());
        otherCopy.Scale(-1.0);
        M_matrix->operator+=(otherCopy);
    }
    return *this;
}

DenseMatrix&
DenseMatrix::
operator*=(const double& coeff)
{
    if (M_matrix)
        M_matrix->Scale(coeff);

    return *this;
}

void
DenseMatrix::
hardCopy(const DenseMatrix& other)
{
    if (other.data())
        M_matrix.reset(new DENSEMATRIX(*other.data()));
}

void
DenseMatrix::
softCopy(const DenseMatrix& other)
{
    M_matrix = other.M_matrix;
}

DenseVector
DenseMatrix::
operator*(const DenseVector& vector)
{
    DenseVector vec;
    if (M_matrix && vector.data())
    {
        std::shared_ptr<DENSEVECTOR> res;
        res.reset(new DENSEVECTOR());
        std::cout << "DenseVector this could give a seg fault" << std::endl << std::flush;
        M_matrix->Multiply(false, *vector.data(), *res);
        vec.data() = res;
    }
    return vec;
}

// void
// DenseMatrix::
// getRowProperty(std::shared_ptr<LifeV::MapEpetra>& outMap)
// {
//     outMap.reset(new LifeV::MapEpetra(M_matrix->domainMap()));
// }
//
// void
// DenseMatrix::
// getColProperty(std::shared_ptr<LifeV::MapEpetra>& outMap)
// {
//     outMap.reset(new LifeV::MapEpetra(M_matrix->rangeMap()));
// }

std::shared_ptr<DENSEMATRIX>&
DenseMatrix::
data()
{
    return M_matrix;
}

std::shared_ptr<DENSEMATRIX>
DenseMatrix::
data() const
{
    return M_matrix;
}

}
