#include "DenseMatrix.hpp"

namespace RedMA
{

DenseMatrix::
DenseMatrix()
{
}

void
DenseMatrix::
add(std::shared_ptr<aMatrix> other)
{
    if (isZero())
    {
       deepCopy(other);
       return;
    }

    if (!other->isZero())
    {
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

std::shared_ptr<aMatrix>
DenseMatrix::
transpose() const
{
    std::shared_ptr<DenseMatrix> retMatrix(new DenseMatrix());

    if (M_matrix)
    {
        unsigned int nrows = M_nCols;
        unsigned int ncols = M_nRows;

        std::shared_ptr<DENSEMATRIX> inMatrix(new DENSEMATRIX());

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
setMatrix(std::shared_ptr<DENSEMATRIX> matrix)
{
    if (matrix)
    {
        M_matrix = matrix;
        this->M_nRows = M_matrix->M();
        this->M_nCols = M_matrix->N();
    }
}

std::shared_ptr<aMatrix>
DenseMatrix::
multiplyByMatrix(std::shared_ptr<aMatrix> other)
{
    std::shared_ptr<DenseMatrix> retMat(new DenseMatrix());

    if (isZero() || other->isZero())
        return retMat;

    std::shared_ptr<DENSEMATRIX> otherMatrix = std::static_pointer_cast<DENSEMATRIX>(other->data());

    if (!isZero() && !other->isZero())
    {
        std::shared_ptr<DENSEMATRIX> mat;
        mat.reset(new DENSEMATRIX(nRows(), other->nCols()));
        mat->Multiply('N', 'N', 1.0, *M_matrix, *otherMatrix, 0.0);
        setMatrix(mat);
    }
}

std::shared_ptr<aVector>
DenseMatrix::
multiplyByVector(std::shared_ptr<aVector> vector)
{
    if (isZero())
    {
        // place holder for zero vector
        return std::shared_ptr<DenseVector>(new DenseVector());
    }

    std::shared_ptr<DENSEVECTOR> otherVector = std::static_pointer_cast<DENSEVECTOR>(vector->data());

    std::shared_ptr<DenseVector> retVec(new DenseVector());
    if (!isZero() && !vector->isZero())
    {
        std::shared_ptr<DENSEVECTOR> res;
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
shallowCopy(std::shared_ptr<aDataWrapper> other)
{
    if (other)
    {
        auto otherMatrix = convert<DenseMatrix>(other);
        std::shared_ptr<DENSEMATRIX> otherMatrixPtr
            (static_cast<DENSEMATRIX*>(otherMatrix->data().get()));
        setMatrix(otherMatrixPtr);
    }
}

void
DenseMatrix::
deepCopy(std::shared_ptr<aDataWrapper> other)
{
    if (other)
    {
        auto otherMatrix = convert<DenseMatrix>(other);
        std::shared_ptr<DENSEMATRIX> otherMatrixPtr
            (static_cast<DENSEMATRIX*>(otherMatrix->data().get()));
        std::shared_ptr<DENSEMATRIX> newMatrix
            (new DENSEMATRIX(*otherMatrixPtr));
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
        std::shared_ptr<DENSEMATRIX> newMatrix
            (new DENSEMATRIX(*M_matrix));
        retMatrix->setMatrix(newMatrix);
    }
    return retMatrix;
}

std::shared_ptr<void>
DenseMatrix::
data() const
{
    return M_matrix;
}

void
DenseMatrix::
setData(std::shared_ptr<void> data)
{
    setData(std::static_pointer_cast<DENSEMATRIX>(data));
}

}
