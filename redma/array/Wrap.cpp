#include "Wrap.hpp"

namespace RedMA
{

shp<DenseMatrix> wrap(shp<DENSEMATRIX> matrix)
{
    shp<DenseMatrix> wrappedMatrix(new DenseMatrix());
    wrappedMatrix->setMatrix(matrix);
    return wrappedMatrix;
}

shp<DenseVector> wrap(shp<DENSEVECTOR> vector)
{
    shp<DenseVector> wrappedVector(new DenseVector());
    wrappedVector->setVector(vector);
    return wrappedVector;
}

shp<SparseMatrix> wrap(shp<MATRIXEPETRA> matrix)
{
    shp<SparseMatrix> wrappedMatrix(new SparseMatrix());
    wrappedMatrix->setMatrix(matrix);
    return wrappedMatrix;
}

shp<DistributedVector> wrap(shp<VECTOREPETRA> vector)
{
    shp<DistributedVector> wrappedVector(new DistributedVector());
    wrappedVector->setVector(vector);
    return wrappedVector;
}

}
