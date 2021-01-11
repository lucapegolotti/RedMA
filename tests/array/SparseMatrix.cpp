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
#include <redma/array/SparseMatrix.hpp>

#include <tests/AtomicTest.hpp>

#include <time.h>

using namespace RedMA;

shp<Epetra_Comm> generateComm()
{
    #ifdef HAVE_MPI
    shp<Epetra_Comm> comm (new Epetra_MpiComm(MPI_COMM_WORLD));
    #else
    shp<Epetra_Comm> comm(new Epetra_SerialComm ());
    #endif
    return comm;
}

shp<LifeV::MapEpetra> generateMap(unsigned int size, shp<Epetra_Comm> comm)
{
    shp<LifeV::MapEpetra> map;
    map.reset(new LifeV::MapEpetra(size, size, 0, comm));

    return map;
}

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
shp<MATRIXEPETRA> fromDense(shp<DENSEMATRIX> otherMatrix)
{
    int N = otherMatrix->M();
    int M = otherMatrix->N();

    shp<Epetra_Comm> comm = generateComm();
    shp<LifeV::MapEpetra> rowMap = generateMap(N, comm);
    shp<LifeV::MapEpetra> colMap = generateMap(M, comm);

    shp<MATRIXEPETRA> mat(new MATRIXEPETRA(*rowMap,M,false));

    for (unsigned int i = 0; i < N; i++)
    {
        for (unsigned int j = 0; j < M; j++)
            mat->addToCoefficient(i, j, otherMatrix->operator()(i,j));
    }
    mat->globalAssemble(colMap, rowMap);
    return mat;
}

shp<MATRIXEPETRA> randmat(unsigned int N, unsigned int M)
{
    return fromDense(randmatdense(N,M));
}

shp<VECTOREPETRA> randvec(unsigned int N)
{
    shp<Epetra_Comm> comm = generateComm();
    shp<LifeV::MapEpetra> map = generateMap(N, comm);
    shp<VECTOREPETRA> vec(new VECTOREPETRA(map));

    for (unsigned int i = 0; i < N; i++)
        vec->operator[](i) = rand() % 10 + 1;

    return vec;
}

bool checkEqual(shp<MATRIXEPETRA> mat1, shp<MATRIXEPETRA> mat2, int numrows, int numcols)
{
    for (unsigned int i = 0; i < numrows; i++)
    {
        int numEntries1;
        int numEntries2;
        double* values1 = new double[numcols];
        double* values2 = new double[numcols];
        int* indices1 = new int[numcols];
        int* indices2 = new int[numcols];
        mat1->matrixPtr()->ExtractGlobalRowCopy(i,numcols,numEntries1,values1,indices1);
        mat2->matrixPtr()->ExtractGlobalRowCopy(i,numcols,numEntries2,values2,indices2);
        if (numEntries1 != numEntries2)
            return false;

        for (unsigned int j = 0; j < numEntries1; j++)
        {
            if (abs(values1[indices1[j]] - values2[indices2[j]]) > TZERO)
                return false;
        }

        delete[] values1;
        delete[] values2;
        delete[] indices1;
        delete[] indices2;
    }
    return true;
}

bool checkEqual(shp<VECTOREPETRA> vec1, shp<VECTOREPETRA> vec2, int L)
{
    for (unsigned int i = 0; i < L; i++)
    {
        if (abs(vec1->operator[](i) - vec2->operator[](i)) > TZERO)
            return false;
    }
    return true;
}

int constructor1()
{
    SparseMatrix dMatrix;

    return SUCCESS;
}

int constructor2()
{
    int N = 10;
    int M = 5;
    std::vector<shp<DistributedVector>> vectors;
    shp<DENSEMATRIX> dmat(new DENSEMATRIX(N,M));

    for (unsigned int i = 0; i < M; i++)
    {
        auto auxvec = randvec(N);
        for (unsigned int j = 0; j < N; j++)
        {
            dmat->operator()(j,i) = auxvec->operator[](j);
        }
        shp<DistributedVector> newVector(new DistributedVector());
        newVector->setVector(auxvec);
        vectors.push_back(newVector);
    }

    shp<SparseMatrix> smat(new SparseMatrix(vectors));
    bool areEqual = checkEqual(fromDense(dmat), smat->getMatrix(), N, M);

    if (areEqual)
        return SUCCESS;

    return FAILURE;
}

int constructor3()
{
    int N = 10;
    int M = 5;
    std::vector<shp<VECTOREPETRA>> vectors;
    shp<DENSEMATRIX> dmat(new DENSEMATRIX(N,M));

    for (unsigned int i = 0; i < M; i++)
    {
        auto auxvec = randvec(N);
        for (unsigned int j = 0; j < N; j++)
        {
            dmat->operator()(j,i) = auxvec->operator[](j);
        }
        vectors.push_back(auxvec);
    }

    shp<SparseMatrix> smat(new SparseMatrix(vectors));
    bool areEqual = checkEqual(fromDense(dmat), smat->getMatrix(), N, M);

    if (areEqual)
        return SUCCESS;

    return FAILURE;
}

int setMatrix()
{
    shp<MATRIXEPETRA> randm = randmat(3,4);

    SparseMatrix dMatrix;
    dMatrix.setMatrix(randm);

    bool areEqual = checkEqual(dMatrix.getMatrix(), randm, 3, 4);

    return !areEqual;
}

// case in which both matrices are filled
int add1()
{
    unsigned int N = 5;
    unsigned int M = 6;

    shp<MATRIXEPETRA> mat1 = randmat(N,M);
    shp<MATRIXEPETRA> mat1Copy(new MATRIXEPETRA(*mat1));
    shp<MATRIXEPETRA> mat2 = randmat(N,M);

    shp<SparseMatrix> dmat1(new SparseMatrix());
    dmat1->setMatrix(mat1);

    shp<SparseMatrix> dmat2(new SparseMatrix());
    dmat2->setMatrix(mat2);

    dmat1->add(dmat2);

    *mat1Copy += *mat2;

    bool areEqual = checkEqual(mat1Copy, dmat1->getMatrix(), N, M);

    if (areEqual)
        return SUCCESS;

    return FAILURE;
}

// non consistent dimensions
int add2()
{
    unsigned int N = 3;
    unsigned int M = 4;

    shp<MATRIXEPETRA> mat1 = randmat(N,M);
    shp<MATRIXEPETRA> mat2 = randmat(N,M+1);

    shp<SparseMatrix> dmat1(new SparseMatrix());
    dmat1->setMatrix(mat1);

    shp<SparseMatrix> dmat2(new SparseMatrix());
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

    shp<SparseMatrix> dmat1(new SparseMatrix());
    shp<SparseMatrix> dmat2(new SparseMatrix());
    shp<MATRIXEPETRA> mat2 = randmat(N,M);

    dmat2->setMatrix(mat2);
    dmat1->add(dmat2);
    areEqual = checkEqual(mat2, dmat1->getMatrix(), N, M);

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

    shp<SparseMatrix> dmat1(new SparseMatrix());
    shp<SparseMatrix> dmat2(new SparseMatrix());
    shp<MATRIXEPETRA> mat1 = randmat(N,M);

    dmat1->setMatrix(mat1);
    dmat1->add(dmat2);

    areEqual = checkEqual(mat1, dmat1->getMatrix(), N, M);

    if (areEqual)
        return SUCCESS;

    return FAILURE;
}

int multiplyByScalar()
{
    unsigned int N = 3;
    unsigned int M = 4;

    shp<SparseMatrix> dmat1(new SparseMatrix());
    shp<MATRIXEPETRA> mat1 = randmat(N,M);

    double normm = mat1->normInf();

    dmat1->setMatrix(mat1);
    dmat1->multiplyByScalar(1.234);

    if (abs(normm * 1.234 - dmat1->getMatrix()->normInf()) < TZERO)
        return SUCCESS;

    return FAILURE;
}

// second matrix is zero
int multiplyByMatrix1()
{
    unsigned int N = 3;
    unsigned int M = 4;

    shp<SparseMatrix> dmat1(new SparseMatrix());
    shp<MATRIXEPETRA> mat1 = randmat(N,M);
    dmat1->setMatrix(mat1);

    shp<SparseMatrix> dmat2(new SparseMatrix());
    shp<SparseMatrix> res = convert<SparseMatrix>(dmat1->multiplyByMatrix(dmat2));

    if (res->isZero())
        return SUCCESS;

    return FAILURE;
}

int multiplyByMatrix2()
{
    unsigned int N = 5;
    unsigned int M = 6;
    unsigned int L = 7;

    shp<SparseMatrix> dmat1(new SparseMatrix());
    shp<DENSEMATRIX> mat1 = randmatdense(N,M);
    dmat1->setMatrix(fromDense(mat1));

    shp<SparseMatrix> dmat2(new SparseMatrix());
    shp<DENSEMATRIX> mat2 = randmatdense(M,L);
    dmat2->setMatrix(fromDense(mat2));

    shp<SparseMatrix> res = convert<SparseMatrix>(dmat1->multiplyByMatrix(dmat2));

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

    bool areEqual = checkEqual(fromDense(thres),res->getMatrix(),N,L);

    if (areEqual)
        return SUCCESS;

    return FAILURE;
}

int multiplyByVector1()
{
    unsigned int N = 3;
    unsigned int M = 5;

    shp<SparseMatrix> dmat(new SparseMatrix());
    shp<MATRIXEPETRA> mat = randmat(N,M);
    dmat->setMatrix(mat);

    shp<DistributedVector> dvec(new DistributedVector());

    shp<DistributedVector> res = convert<DistributedVector>(dmat->multiplyByVector(dvec));

    if (res->isZero())
        return SUCCESS;

    return FAILURE;
}

int multiplyByVector2()
{
    unsigned int N = 3;
    unsigned int M = 5;

    shp<SparseMatrix> dmat(new SparseMatrix());
    shp<DENSEMATRIX> mat = randmatdense(N,M);
    dmat->setMatrix(fromDense(mat));

    shp<DistributedVector> dvec(new DistributedVector());
    shp<VECTOREPETRA> vec = randvec(M);
    dvec->setVector(vec);

    shp<DistributedVector> res = convert<DistributedVector>(dmat->multiplyByVector(dvec));

    shp<Epetra_Comm> comm = generateComm();
    shp<LifeV::MapEpetra> map = generateMap(N, comm);
    shp<VECTOREPETRA> thres(new VECTOREPETRA(map));
    for (unsigned int i = 0; i < N; i++)
    {
        double par = 0;
        for (unsigned int j = 0; j < M; j++)
        {
            par += mat->operator()(i,j) * vec->operator[](j);
        }
        thres->operator[](i) = par;
    }

    bool areEqual = checkEqual(thres, res->getVector(), N);

    if (areEqual)
        return SUCCESS;

    return FAILURE;
}

// we just check that no errors are thrown
int _clone()
{
    unsigned int M = 5;

    shp<SparseMatrix> dmat(new SparseMatrix());
    shp<MATRIXEPETRA> mat = randmat(M,M);
    dmat->setMatrix(mat);

    shp<aMatrix> aux(dmat->clone());

    return SUCCESS;
}

int transpose()
{
    unsigned int N = 6;
    unsigned int M = 5;

    shp<SparseMatrix> dmat(new SparseMatrix());
    shp<DENSEMATRIX> mat = randmatdense(N,M);
    dmat->setMatrix(fromDense(mat));

    shp<SparseMatrix> dmatT = convert<SparseMatrix>(dmat->transpose());
    shp<MATRIXEPETRA> matT = dmatT->getMatrix();

    shp<DENSEMATRIX> densematT(new DENSEMATRIX(M,N));
    for (unsigned int i = 0; i < N; i++)
    {
        for (unsigned int j = 0; j < M; j++)
        {
            densematT->operator()(j,i) = mat->operator()(i,j);
        }
    }

    bool areEqual = checkEqual(fromDense(densematT), matT, M, N);

    if (areEqual)
        return SUCCESS;

    return FAILURE;
}

int shallowCopy()
{
    unsigned int N = 6;
    unsigned int M = 5;

    shp<SparseMatrix> dmat(new SparseMatrix());
    dmat->setMatrix(randmat(N,M));

    shp<SparseMatrix> otherdmat(new SparseMatrix());
    otherdmat->shallowCopy(dmat);

    otherdmat->multiplyByScalar(2);

    // check if also dmat was modified
    if (otherdmat->getMatrix()->normInf() == dmat->getMatrix()->normInf())
        return SUCCESS;

    return FAILURE;
}

int deepCopy()
{
    unsigned int N = 6;
    unsigned int M = 5;

    shp<SparseMatrix> dmat(new SparseMatrix());
    dmat->setMatrix(randmat(N,M));

    shp<SparseMatrix> otherdmat(new SparseMatrix());
    otherdmat->deepCopy(dmat);

    otherdmat->multiplyByScalar(2);

    // check if dmat wasn't modified
    if (otherdmat->getMatrix()->normInf() != dmat->getMatrix()->normInf())
        return SUCCESS;

    return FAILURE;
}

int toDenseMatrix()
{
    unsigned int N = 6;
    unsigned int M = 5;

    shp<SparseMatrix> smat(new SparseMatrix());
    shp<DENSEMATRIX> mat = randmatdense(N,M);
    smat->setMatrix(fromDense(mat));

    auto dmat = smat->toDenseMatrix();
    dmat.multiplyByScalar(-1);
    *dmat.getMatrix() += *mat;
    double normDiff = dmat.getMatrix()->NormInf();

    if (normDiff < TZERO)
        return SUCCESS;

    return FAILURE;
}

int convertDenseMatrix()
{
    int N = 5;
    int M = 6;

    auto comm = generateComm();
    auto dmat = randmatdense(N,M);
    shp<DenseMatrix> dmat2(new DenseMatrix());
    dmat2->setMatrix(dmat);

    auto smat = SparseMatrix::convertDenseMatrix(dmat2,comm);

    bool areEqual = checkEqual(smat->getMatrix(), fromDense(dmat), N, M);

    if (areEqual)
        return SUCCESS;

    return FAILURE;
}

int main()
{
    MPI_Init(nullptr, nullptr);
    srand(time(NULL));

    int status = SUCCESS;

    status |= AtomicTest("Constructor1", &constructor1).run();
    status |= AtomicTest("Constructor2", &constructor2).run();
    status |= AtomicTest("Constructor3", &constructor3).run();
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
    status |= AtomicTest("ToDenseMatrix", &toDenseMatrix).run();
    status |= AtomicTest("ConvertDenseMatrix", &convertDenseMatrix).run();

    return status;
}
