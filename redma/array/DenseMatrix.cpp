#include "DenseMatrix.hpp"

namespace RedMA
{

DenseMatrix::
DenseMatrix()
{
}

void
DenseMatrix::
add(shp<aMatrix> other)
{
    if (isZero())
    {
       deepCopy(other);
       return;
    }

    if (!other->isZero())
    {
        if (other->nRows() != nRows() || other->nCols() != nCols())
            throw new Exception("[DenseMatrix::add] inconsistent dimensions of matrices");

        shp<DenseMatrix> otherMatrix = convert<DenseMatrix>(other);
        (*M_matrix) += *static_cast<DENSEMATRIX*>(otherMatrix->data().get());
    }
}

void
DenseMatrix::
multiplyByScalar(const double& coeff)
{
    if (!isZero())
        M_matrix->Scale(coeff);
}

shp<aMatrix>
DenseMatrix::
transpose() const
{
    shp<DenseMatrix> retMatrix(new DenseMatrix());

    if (M_matrix)
    {
        unsigned int nrows = M_nCols;
        unsigned int ncols = M_nRows;

        shp<DENSEMATRIX> inMatrix(new DENSEMATRIX());

        inMatrix->Reshape(nrows,ncols);

        for (unsigned int i = 0; i < nrows; i++)
        {
            for (unsigned int j = 0; j < ncols; j++)
            {
                (*inMatrix)(i,j) = (*M_matrix)(j,i);
            }
        }

        retMatrix->setMatrix(inMatrix);
    }

    return retMatrix;
}

void
DenseMatrix::
setMatrix(shp<DENSEMATRIX> matrix)
{
    if (matrix)
    {
        M_matrix = matrix;
        this->M_nRows = M_matrix->M();
        this->M_nCols = M_matrix->N();
    }
}

shp<aMatrix>
DenseMatrix::
multiplyByMatrix(shp<aMatrix> other)
{
    shp<DenseMatrix> retMat(new DenseMatrix());

    if (isZero() || other->isZero())
        return retMat;

    shp<DENSEMATRIX> otherMatrix = spcast<DENSEMATRIX>(other->data());

    if (!isZero() && !other->isZero())
    {
        shp<DENSEMATRIX> mat;
        mat.reset(new DENSEMATRIX(nRows(), other->nCols()));
        mat->Multiply('N', 'N', 1.0, *M_matrix, *otherMatrix, 0.0);
        retMat->setMatrix(mat);
    }

    return retMat;
}

shp<aVector>
DenseMatrix::
multiplyByVector(shp<aVector> vector)
{
    shp<DenseVector> retVec(new DenseVector());

    if (isZero())
        return retVec;

    shp<DENSEVECTOR> otherVector = spcast<DENSEVECTOR>(vector->data());

    if (!isZero() && !vector->isZero())
    {
        if (vector->nRows() != nCols())
            throw new Exception("[DenseMatrix::multiplyByVector] inconsistent dimensions!");
        shp<DENSEVECTOR> res;
        res.reset(new DENSEVECTOR(M_nRows));
        M_matrix->Multiply(false, *otherVector, *res);
        retVec->setVector(res);
    }

    return retVec;
}

void
DenseMatrix::
dump(std::string filename) const
{
    std::ofstream outfile(filename);

    if (M_matrix)
    {
        unsigned int M = M_matrix->M();
        unsigned int N = M_matrix->N();


        for (unsigned int i = 0; i < M; i++)
        {
            for (unsigned int j = 0; j < N; j++)
            {
                outfile << (*M_matrix)(i,j);
                if (j != N-1)
                    outfile << ",";
            }
            outfile << "\n";
        }
    }

    outfile.close();
}

bool
DenseMatrix::
isZero() const
{
    return M_matrix == nullptr;
}

void
DenseMatrix::
shallowCopy(shp<aDataWrapper> other)
{
    if (other)
    {
        auto otherMatrix = convert<DenseMatrix>(other);
        setMatrix(otherMatrix->getMatrix());
    }
}

void
DenseMatrix::
deepCopy(shp<aDataWrapper> other)
{
    if (other)
    {
        auto otherMatrix = convert<DenseMatrix>(other);
        shp<DENSEMATRIX> otherMatrixPtr = otherMatrix->getMatrix();
        shp<DENSEMATRIX> newMatrix(new DENSEMATRIX(*otherMatrixPtr));
        setMatrix(newMatrix);
    }
}

DenseMatrix*
DenseMatrix::
clone() const
{
    DenseMatrix* retMatrix = new DenseMatrix();
    if (M_matrix)
    {
        shp<DENSEMATRIX> newMatrix
            (new DENSEMATRIX(*M_matrix));
        retMatrix->setMatrix(newMatrix);
    }
    return retMatrix;
}

shp<void>
DenseMatrix::
data() const
{
    return M_matrix;
}

void
DenseMatrix::
setData(shp<void> data)
{
    spcast<DENSEMATRIX>(data);
}

shp<DENSEMATRIX>
DenseMatrix::
getMatrix()
{
    return M_matrix;
}

}
