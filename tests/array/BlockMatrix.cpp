// Reduced Modeling of Arteries (RedMA)
// Copyright (C) 2019  Luca Pegolotti
//
// RedMA is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// RedMA is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with this program.  If not, see <http://www.gnu.org/licenses/>.

#include <redma/RedMA.hpp>
#include <redma/array/BlockMatrix.hpp>

#include <tests/AtomicTest.hpp>

#include <time.h>

using namespace RedMA;

// we test only block matrices composed of dense structures for simplicity
shp<DENSEMATRIX> randmatdense(unsigned int N, unsigned int M)
{
    shp<DENSEMATRIX> mat(new DENSEMATRIX(N,M));

    for (unsigned int i = 0; i < N; i++)
    {
        for (unsigned int j = 0; j < M; j++)
        {
            mat->operator()(i,j) = rand() % 10;
        }
    }
    return mat;
}

shp<DENSEVECTOR> randvecdense(unsigned int N)
{
    shp<DENSEVECTOR> vec(new DENSEVECTOR(N));

    for (unsigned int i = 0; i < N; i++)
        vec->operator()(i) = rand() % 10;

    return vec;
}

shp<BlockMatrix> generateBlockMatrix(unsigned int N, unsigned int M,
                                     std::vector<int> sizesrows,
                                     std::vector<int> sizescols)
{
    shp<BlockMatrix> retMat(new BlockMatrix(N,M));

    for (unsigned int i = 0; i < N; i++)
    {
        for (unsigned int j = 0; j < M; j++)
        {
            retMat->setBlock(i,j,wrap(randmatdense(sizesrows[i],sizescols[j])));
        }
    }

    return retMat;
}

shp<BlockVector> generateBlockVector(unsigned int N, std::vector<int> sizesrows)
{
    shp<BlockVector> retVec(new BlockVector(N));

    for (unsigned int i = 0; i < N; i++)
        retVec->setBlock(i,wrap(randvecdense(sizesrows[i])));

    return retVec;
}

void fillrowscols(int N, int M, std::vector<int>& rows, std::vector<int>& cols)
{
    for (unsigned int i = 0; i < N; i++)
        rows.push_back(i+1);

    for (unsigned int i = M; i > 0; i--)
        cols.push_back(i);
}

int constructor1()
{
    BlockMatrix bMatrix;

    return SUCCESS;
}

int constructor2()
{
    BlockMatrix bMatrix(3,3);

    if (bMatrix.nRows() != 3 || bMatrix.nCols() != 3)
        return FAILURE;

    return SUCCESS;
}

// one dimension is wrong
int add1()
{
    unsigned int N = 5;
    unsigned int M = 6;

    std::vector<int> rows;
    std::vector<int> cols;
    fillrowscols(N,M,rows,cols);

    shp<BlockMatrix> bmat1 = generateBlockMatrix(N, M, rows, cols);
    cols[3] = 1;
    shp<BlockMatrix> bmat2 = generateBlockMatrix(N, M, rows, cols);

    try
    {
        bmat1->add(bmat2);
    }
    catch (Exception* e)
    {
        return SUCCESS;
    }

    return FAILURE;
}

int add2()
{
    unsigned int N = 5;
    unsigned int M = 6;

    std::vector<int> rows;
    std::vector<int> cols;
    fillrowscols(N,M,rows,cols);

    shp<BlockMatrix> bmat1 = generateBlockMatrix(N, M, rows, cols);
    shp<BlockMatrix> bmat2 = generateBlockMatrix(N, M, rows, cols);

    shp<DenseMatrix> block(new DenseMatrix());
    block->deepCopy(bmat1->block(1,1));
    block->add(bmat2->block(1,1));

    bmat1->add(bmat2);

    if (abs(block->getMatrix()->NormInf() - convert<DenseMatrix>(bmat1->block(1,1))->getMatrix()->NormInf()) < TZERO)
        return SUCCESS;

    return FAILURE;
}

int add3()
{
    unsigned int N = 5;
    unsigned int M = 6;

    std::vector<int> rows;
    std::vector<int> cols;
    fillrowscols(N,M,rows,cols);

    shp<BlockMatrix> bmat1(new BlockMatrix());
    shp<BlockMatrix> bmat2 = generateBlockMatrix(N, M, rows, cols);

    bmat1->add(bmat2);

    if (abs(convert<DenseMatrix>(bmat1->block(1,1))->getMatrix()->NormInf() -
            convert<DenseMatrix>(bmat2->block(1,1))->getMatrix()->NormInf()) < TZERO)
        return SUCCESS;

    return FAILURE;
}

// other matrix is empty
int add4()
{
    unsigned int N = 5;
    unsigned int M = 6;

    std::vector<int> rows;
    std::vector<int> cols;
    fillrowscols(N,M,rows,cols);

    shp<BlockMatrix> bmat1 = generateBlockMatrix(N, M, rows, cols);
    shp<BlockMatrix> bmat2(new BlockMatrix());

    double nnorm = convert<DenseMatrix>(bmat1->block(1,1))->getMatrix()->NormInf();
    bmat1->add(bmat2);

    if (abs(nnorm - convert<DenseMatrix>(bmat1->block(1,1))->getMatrix()->NormInf()) < TZERO)
        return SUCCESS;

    return FAILURE;
}

int multiplyByScalar()
{
    unsigned int N = 5;
    unsigned int M = 6;

    std::vector<int> rows;
    std::vector<int> cols;
    fillrowscols(N,M,rows,cols);

    shp<BlockMatrix> bmat = generateBlockMatrix(N, M, rows, cols);

    double nnorm = convert<DenseMatrix>(bmat->block(1,1))->getMatrix()->NormInf();

    bmat->multiplyByScalar(1.234);

    if (abs(nnorm * 1.234 - convert<DenseMatrix>(bmat->block(1,1))->getMatrix()->NormInf()) < TZERO)
        return SUCCESS;

    return FAILURE;
}

int multiplyByMatrix()
{
    unsigned int N = 3;
    unsigned int M = 4;
    unsigned int L = 5;

    std::vector<int> rows1;
    std::vector<int> cols1;
    std::vector<int> rows2;
    std::vector<int> cols2;

    for (unsigned int i = 0; i < N; i++)
        rows1.push_back(2);

    for (unsigned int i = 0; i < M; i++)
        cols1.push_back(2);

    for (unsigned int i = 0; i < M; i++)
        rows2.push_back(2);

    for (unsigned int i = 0; i < L; i++)
        cols2.push_back(2);

    shp<BlockMatrix> bmat1 = generateBlockMatrix(N,M,rows1,cols1);
    shp<BlockMatrix> bmat2 = generateBlockMatrix(M,L,rows2,cols2);

    bmat1->multiplyByMatrix(bmat2);

    return SUCCESS;
}

int multiplyByVector1()
{
    unsigned int N = 3;
    unsigned int M = 5;
    std::vector<int> rows;
    std::vector<int> cols;
    fillrowscols(N,M,rows,cols);

    shp<BlockMatrix> bmat = generateBlockMatrix(N,M,rows,cols);

    shp<BlockVector> bvec(new BlockVector());

    shp<BlockVector> res = convert<BlockVector>(bmat->multiplyByVector(bvec));

    if (res->isZero())
        return SUCCESS;

    return FAILURE;
}

int multiplyByVector2()
{
    unsigned int N = 10;
    unsigned int M = 9;
    std::vector<int> rows;
    std::vector<int> cols;
    fillrowscols(N,M,rows,cols);

    shp<BlockMatrix> bmat = generateBlockMatrix(N,M,rows,cols);
    shp<BlockVector> bvec = generateBlockVector(M,cols);
    shp<BlockVector> res = convert<BlockVector>(bmat->multiplyByVector(bvec));

    // we check only one block
    shp<DenseVector> dvec = convert<DenseVector>(bmat->block(2,0)->multiplyByVector(bvec->block(0)));

    for (unsigned int i = 1; i < M; i++)
        dvec->add(bmat->block(2,i)->multiplyByVector(bvec->block(i)));

    if (abs(convert<DenseVector>(res->block(2))->getVector()->NormInf() - dvec->getVector()->NormInf()) < TZERO)
        return SUCCESS;

    return FAILURE;
}

int multiplyByVector3()
{
    unsigned int N = 3;
    unsigned int M = 5;
    std::vector<int> rows;
    std::vector<int> cols;
    fillrowscols(N,M,rows,cols);

    shp<BlockMatrix> bmat = generateBlockMatrix(N,M,rows,cols);
    cols[1] = 10;
    shp<BlockVector> bvec = generateBlockVector(M,cols);

    try
    {
        shp<BlockVector> res = convert<BlockVector>(bmat->multiplyByVector(bvec));
    }
    catch (Exception* e)
    {
        return SUCCESS;
    }

    return FAILURE;
}


// we just check that no errors are thrown
int clone()
{
    unsigned int N = 3;
    unsigned int M = 5;
    std::vector<int> rows;
    std::vector<int> cols;
    fillrowscols(N,M,rows,cols);

    shp<BlockMatrix> bmat = generateBlockMatrix(N,M,rows,cols);

    shp<aMatrix> aux(bmat->clone());

    return SUCCESS;
}

int transpose()
{
    unsigned int N = 3;
    unsigned int M = 5;
    std::vector<int> rows;
    std::vector<int> cols;
    fillrowscols(N,M,rows,cols);

    shp<BlockMatrix> bmat = generateBlockMatrix(N,M,rows,cols);

    bmat->transpose();

    return SUCCESS;
}

int shallowCopy()
{
    unsigned int N = 3;
    unsigned int M = 5;
    std::vector<int> rows;
    std::vector<int> cols;
    fillrowscols(N,M,rows,cols);

    shp<BlockMatrix> bmat = generateBlockMatrix(N,M,rows,cols);

    shp<BlockMatrix> otherbmat(new BlockMatrix());
    otherbmat->shallowCopy(bmat);

    shp<DENSEMATRIX> m1 = convert<DenseMatrix>(otherbmat->block(2,2))->getMatrix();
    shp<DENSEMATRIX> m2 = convert<DenseMatrix>(bmat->block(2,2))->getMatrix();

    m1->operator()(1,1) = m1->operator()(1,1) + 1;

    // check if also dmat was modified
    if (m1->operator()(1,1) == m2->operator()(1,1))
        return SUCCESS;

    return FAILURE;
}

int deepCopy()
{
    unsigned int N = 3;
    unsigned int M = 5;
    std::vector<int> rows;
    std::vector<int> cols;
    fillrowscols(N,M,rows,cols);

    shp<BlockMatrix> bmat = generateBlockMatrix(N,M,rows,cols);

    shp<BlockMatrix> otherbmat(new BlockMatrix());
    otherbmat->deepCopy(bmat);

    shp<DENSEMATRIX> m1 = convert<DenseMatrix>(otherbmat->block(2,2))->getMatrix();
    shp<DENSEMATRIX> m2 = convert<DenseMatrix>(bmat->block(2,2))->getMatrix();

    m1->operator()(1,1) = m1->operator()(1,1) + 1;

    // check if also dmat was modified
    if (m1->operator()(1,1) != m2->operator()(1,1))
        return SUCCESS;

    return FAILURE;
}

int getSubmatrix()
{
    unsigned int N = 5;
    unsigned int M = 6;
    std::vector<int> rows;
    std::vector<int> cols;
    fillrowscols(N,M,rows,cols);
    shp<BlockMatrix> bmat = generateBlockMatrix(N,M,rows,cols);

    auto bmatsub = bmat->getSubmatrix(2,4,2,4);

    if (bmatsub->nRows() != 3 || bmatsub->nCols() != 3)
        return FAILURE;

    auto b1 = convert<DenseMatrix>(bmatsub->block(0,0))->getMatrix();
    auto b2 = convert<DenseMatrix>(bmat->block(2,2))->getMatrix();

    if (abs(b1->NormInf() - b2->NormInf()) > TZERO)
        return FAILURE;

    return SUCCESS;
}

int main()
{
    srand(time(NULL));

    int status = SUCCESS;

    status |= AtomicTest("Constructor1", &constructor1).run();
    status |= AtomicTest("Constructor2", &constructor2).run();
    status |= AtomicTest("Add1", &add1).run();
    status |= AtomicTest("Add2", &add2).run();
    status |= AtomicTest("Add3", &add3).run();
    status |= AtomicTest("Add4", &add4).run();
    status |= AtomicTest("MultiplyByScalar", &multiplyByScalar).run();
    status |= AtomicTest("MultiplyByMatrix", &multiplyByMatrix).run();
    status |= AtomicTest("MultiplyByVector1", &multiplyByVector1).run();
    status |= AtomicTest("MultiplyByVector2", &multiplyByVector2).run();
    status |= AtomicTest("MultiplyByVector3", &multiplyByVector3).run();
    status |= AtomicTest("Clone", &clone).run();
    status |= AtomicTest("Transpose", &transpose).run();
    status |= AtomicTest("ShallowCopy", &shallowCopy).run();
    status |= AtomicTest("DeepCopy", &deepCopy).run();
    status |= AtomicTest("Submatrix", &getSubmatrix).run();

    return status;
}
