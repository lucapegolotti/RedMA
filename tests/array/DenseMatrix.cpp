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
#include <redma/array/DenseMatrix.hpp>

#include <tests/AtomicTest.hpp>

#include <time.h>

using namespace RedMA;

shp<DENSEMATRIX> identity(unsigned int N)
{
    shp<DENSEMATRIX> mat(new DENSEMATRIX(N,N));

    for (unsigned int i = 0; i < N; i++)
    {
        for (unsigned int j = 0; j < N; j++)
        {
            mat->operator()(i,j) = (i == j);
        }
    }
    return mat;
}

shp<DENSEMATRIX> randmat(unsigned int N, unsigned int M)
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

shp<DENSEVECTOR> randvec(unsigned int N)
{
    shp<DENSEVECTOR> vec(new DENSEVECTOR(N));

    for (unsigned int i = 0; i < N; i++)
        vec->operator()(i) = rand() % 10;

    return vec;
}

bool checkEqual(shp<DENSEMATRIX> mat1, shp<DENSEMATRIX> mat2)
{
    for (unsigned int i = 0; i < mat1->M(); i++)
    {
        for (unsigned int j = 0; j < mat1->N(); j++)
        {
            if (abs(mat1->operator()(i,j) - mat2->operator()(i,j)) > TZERO)
                return false;
        }
    }
    return true;
}

bool checkEqual(shp<DENSEVECTOR> vec1, shp<DENSEVECTOR> vec2)
{
    for (unsigned int i = 0; i < vec1->Length(); i++)
    {
        if (abs(vec1->operator()(i) - vec2->operator()(i)) > TZERO)
            return false;
    }
    return true;
}

int constructor()
{
    DenseMatrix dMatrix;

    return SUCCESS;
}

int setMatrix()
{
    DenseMatrix dMatrix;
    dMatrix.setMatrix(identity(3));

    bool areEqual = checkEqual(dMatrix.getMatrix(), identity(3));

    return !areEqual;
}

// case in which both matrices are filled
int add1()
{
    unsigned int N = 3;
    unsigned int M = 4;

    shp<DENSEMATRIX> mat1 = randmat(N,M);
    shp<DENSEMATRIX> mat1Copy(new DENSEMATRIX(*mat1));
    shp<DENSEMATRIX> mat2 = randmat(N,M);

    shp<DenseMatrix> dmat1(new DenseMatrix());
    dmat1->setMatrix(mat1);

    shp<DenseMatrix> dmat2(new DenseMatrix());
    dmat2->setMatrix(mat2);

    dmat1->add(dmat2);

    for (unsigned int i = 0; i < N; i++)
    {
        for (unsigned int j = 0; j < M; j++)
        {
            double expected = mat1Copy->operator()(i,j) + mat2->operator()(i,j);
            if (abs(dmat1->getMatrix()->operator()(i,j) - expected) > TZERO)
                return FAILURE;
        }
    }
    return SUCCESS;
}

// non consistent dimensions
int add2()
{
    unsigned int N = 3;
    unsigned int M = 4;

    shp<DENSEMATRIX> mat1 = randmat(N,M);
    shp<DENSEMATRIX> mat2 = randmat(N,M+1);

    shp<DenseMatrix> dmat1(new DenseMatrix());
    dmat1->setMatrix(mat1);

    shp<DenseMatrix> dmat2(new DenseMatrix());
    dmat2->setMatrix(mat2);

    try
    {
        dmat1->add(dmat2);
    }
    catch (Exception* e)
    {
        return SUCCESS;
    }

    return FAILURE;
}

// current matrix is empty
int add3()
{
    bool areEqual;

    unsigned int N = 3;
    unsigned int M = 4;

    shp<DenseMatrix> dmat1(new DenseMatrix());
    shp<DenseMatrix> dmat2(new DenseMatrix());
    shp<DENSEMATRIX> mat2 = randmat(N,M);

    dmat2->setMatrix(mat2);
    dmat1->add(dmat2);
    areEqual = checkEqual(mat2, dmat1->getMatrix());

    if (areEqual)
        return SUCCESS;

    return FAILURE;
}

// other matrix is empty
int add4()
{
    bool areEqual;

    unsigned int N = 3;
    unsigned int M = 4;

    shp<DenseMatrix> dmat1(new DenseMatrix());
    shp<DenseMatrix> dmat2(new DenseMatrix());
    shp<DENSEMATRIX> mat1 = randmat(N,M);

    dmat1->setMatrix(mat1);
    dmat1->add(dmat2);

    areEqual = checkEqual(mat1, dmat1->getMatrix());

    if (areEqual)
        return SUCCESS;

    return FAILURE;
}

int multiplyByScalar()
{
    unsigned int N = 3;
    unsigned int M = 4;

    shp<DenseMatrix> dmat1(new DenseMatrix());
    shp<DENSEMATRIX> mat1 = randmat(N,M);

    double normm = mat1->NormInf();

    dmat1->setMatrix(mat1);
    dmat1->multiplyByScalar(1.234);

    if (abs(normm * 1.234 - dmat1->getMatrix()->NormInf()) < TZERO)
        return SUCCESS;

    return FAILURE;
}

// second matrix is zero
int multiplyByMatrix1()
{
    unsigned int N = 3;
    unsigned int M = 4;

    shp<DenseMatrix> dmat1(new DenseMatrix());
    shp<DENSEMATRIX> mat1 = randmat(N,M);
    dmat1->setMatrix(mat1);

    shp<DenseMatrix> dmat2(new DenseMatrix());
    shp<DenseMatrix> res = convert<DenseMatrix>(dmat1->multiplyByMatrix(dmat2));

    if (res->isZero())
        return SUCCESS;

    return FAILURE;
}

int multiplyByMatrix2()
{
    unsigned int N = 3;
    unsigned int M = 4;
    unsigned int L = 5;

    shp<DenseMatrix> dmat1(new DenseMatrix());
    shp<DENSEMATRIX> mat1 = randmat(N,M);
    dmat1->setMatrix(mat1);

    shp<DenseMatrix> dmat2(new DenseMatrix());
    shp<DENSEMATRIX> mat2 = randmat(M,L);
    dmat2->setMatrix(mat2);

    shp<DenseMatrix> res = convert<DenseMatrix>(dmat1->multiplyByMatrix(dmat2));

    shp<DENSEMATRIX> thres(new DENSEMATRIX(N,L));
    for (unsigned int i = 0; i < N; i++)
    {
        for (unsigned int j = 0; j < L; j++)
        {
            double par = 0;
            for (unsigned int k = 0; k < M; k++)
            {
                par += mat1->operator()(i,k) * mat2->operator()(k,j);
            }
            thres->operator()(i,j) = par;
        }
    }

    bool areEqual = checkEqual(thres, res->getMatrix());

    if (areEqual)
        return SUCCESS;

    return FAILURE;
}

int multiplyByVector1()
{
    unsigned int N = 3;
    unsigned int M = 5;

    shp<DenseMatrix> dmat(new DenseMatrix());
    shp<DENSEMATRIX> mat = randmat(N,M);
    dmat->setMatrix(mat);

    shp<DenseVector> dvec(new DenseVector());

    shp<DenseVector> res = convert<DenseVector>(dmat->multiplyByVector(dvec));

    if (res->isZero())
        return SUCCESS;

    return FAILURE;
}

int multiplyByVector2()
{
    unsigned int N = 3;
    unsigned int M = 5;

    shp<DenseMatrix> dmat(new DenseMatrix());
    shp<DENSEMATRIX> mat = randmat(N,M);
    dmat->setMatrix(mat);

    shp<DenseVector> dvec(new DenseVector());
    shp<DENSEVECTOR> vec = randvec(M);
    dvec->setVector(vec);

    shp<DenseVector> res = convert<DenseVector>(dmat->multiplyByVector(dvec));

    shp<DENSEVECTOR> thres(new DENSEVECTOR(N));
    for (unsigned int i = 0; i < N; i++)
    {
        double par = 0;
        for (unsigned int j = 0; j < M; j++)
        {
            par += mat->operator()(i,j) * vec->operator()(j);
        }
        thres->operator()(i) = par;
    }

    bool areEqual = checkEqual(thres, res->getVector());

    if (areEqual)
        return SUCCESS;

    return FAILURE;
}

// we just check that no errors are thrown
int _clone()
{
    unsigned int M = 5;

    shp<DenseMatrix> dmat(new DenseMatrix());
    shp<DENSEMATRIX> mat = randmat(M,M);
    dmat->setMatrix(mat);

    shp<aMatrix> aux(dmat->clone());

    return SUCCESS;
}

int transpose()
{
    unsigned int N = 6;
    unsigned int M = 5;

    shp<DenseMatrix> dmat(new DenseMatrix());
    shp<DENSEMATRIX> mat = randmat(N,M);
    dmat->setMatrix(mat);

    shp<DenseMatrix> dmatT = convert<DenseMatrix>(dmat->transpose());
    shp<DENSEMATRIX> matT = dmatT->getMatrix();

    for (unsigned int i = 0; i < N; i++)
    {
        for (unsigned int j = 0; j < M; j++)
        {
            if (abs((*mat)(i,j)-(*matT)(j,i)) > TZERO)
                return FAILURE;
        }
    }
    return SUCCESS;
}

int shallowCopy()
{
    unsigned int N = 6;
    unsigned int M = 5;

    shp<DenseMatrix> dmat(new DenseMatrix());
    dmat->setMatrix(randmat(N,M));

    shp<DenseMatrix> otherdmat(new DenseMatrix());
    otherdmat->shallowCopy(dmat);

    otherdmat->getMatrix()->operator()(1,1) = otherdmat->getMatrix()->operator()(1,1) + 1;

    // check if also dmat was modified
    if (otherdmat->getMatrix()->operator()(1,1) == dmat->getMatrix()->operator()(1,1))
        return SUCCESS;

    return FAILURE;
}

int deepCopy()
{
    unsigned int N = 6;
    unsigned int M = 5;

    shp<DenseMatrix> dmat(new DenseMatrix());
    dmat->setMatrix(randmat(N,M));

    shp<DenseMatrix> otherdmat(new DenseMatrix());
    otherdmat->deepCopy(dmat);

    otherdmat->getMatrix()->operator()(1,1) = otherdmat->getMatrix()->operator()(1,1) + 1;

    // check if dmat wasn't modified
    if (otherdmat->getMatrix()->operator()(1,1) != dmat->getMatrix()->operator()(1,1))
        return SUCCESS;

    return FAILURE;
}

int main()
{
    srand(time(NULL));

    int status = SUCCESS;

    status |= AtomicTest("Constructor", &constructor).run();
    status |= AtomicTest("SetMatrix", &setMatrix).run();
    status |= AtomicTest("Add1", &add1).run();
    status |= AtomicTest("Add2", &add2).run();
    status |= AtomicTest("Add3", &add3).run();
    status |= AtomicTest("Add4", &add4).run();
    status |= AtomicTest("MultiplyByScalar", &multiplyByScalar).run();
    status |= AtomicTest("MultiplyByMatrix1", &multiplyByMatrix1).run();
    status |= AtomicTest("MultiplyByMatrix2", &multiplyByMatrix2).run();
    status |= AtomicTest("MultiplyByVector1", &multiplyByVector1).run();
    status |= AtomicTest("MultiplyByVector2", &multiplyByVector1).run();
    status |= AtomicTest("Clone", &_clone).run();
    status |= AtomicTest("Transpose", &transpose).run();
    status |= AtomicTest("ShallowCopy", &shallowCopy).run();
    status |= AtomicTest("DeepCopy", &deepCopy).run();

    return status;
}
