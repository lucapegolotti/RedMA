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
#include <redma/array/BlockVector.hpp>

#include <tests/AtomicTest.hpp>

#include <time.h>

using namespace RedMA;

shp<Epetra_Comm> generateComm()
{
    #ifdef HAVE_MPI
    std::shared_ptr<Epetra_Comm> comm (new Epetra_MpiComm(MPI_COMM_WORLD));
    #else
    std::shared_ptr<Epetra_Comm> comm(new Epetra_SerialComm ());
    #endif
    return comm;
}

shp<LifeV::MapEpetra> generateMap(unsigned int size, shp<Epetra_Comm> comm)
{
    shp<LifeV::MapEpetra> map;
    map.reset(new LifeV::MapEpetra(size, size, 0, comm));

    return map;
}

shp<VECTOREPETRA> randvecepetra(unsigned int N)
{
    shp<Epetra_Comm> comm = generateComm();
    shp<LifeV::MapEpetra> map = generateMap(N, comm);
    shp<VECTOREPETRA> vec(new VECTOREPETRA(map));

    for (unsigned int i = 0; i < N; i++)
        vec->operator[](i) = rand() % 10;

    return vec;
}

shp<DENSEVECTOR> randvecdense(unsigned int N)
{
    shp<DENSEVECTOR> vec(new DENSEVECTOR(N));

    for (unsigned int i = 0; i < N; i++)
        vec->operator()(i) = rand() % 10;

    return vec;
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
    BlockVector bVector;

    return SUCCESS;
}

int constructor2()
{
    BlockVector bVector(3);

    if (bVector.nRows() != 3)
        return FAILURE;

    return SUCCESS;
}

// case in which both vectors are filled. We just check the norm of the result for
// simplicity
int add1()
{
    unsigned int N = 4;
    //
    // shp<VECTOREPETRA> vec1 = randvec(N);
    // shp<VECTOREPETRA> vec1Copy(new VECTOREPETRA(*vec1));
    // shp<VECTOREPETRA> vec2 = randvec(N);

    shp<BlockVector> bvec1(new BlockVector(2));
    shp<VECTOREPETRA> vec11 = randvecepetra(N);
    shp<DENSEVECTOR> vec12 = randvecdense(N);
    bvec1->setBlock(0, wrap(vec11));
    bvec1->setBlock(1, wrap(vec12));

    shp<VECTOREPETRA> vec11copy(new VECTOREPETRA(*vec11));
    shp<DENSEVECTOR> vec12copy(new DENSEVECTOR(*vec12));

    shp<BlockVector> bvec2(new BlockVector(2));
    shp<VECTOREPETRA> vec21 = randvecepetra(N);
    shp<DENSEVECTOR> vec22 = randvecdense(N);
    bvec2->setBlock(0, wrap(vec21));
    bvec2->setBlock(1, wrap(vec22));

    *vec11copy += *vec21;
    *vec12copy += *vec22;

    double expectedNorm = 0;
    expectedNorm += vec11copy->norm2() * vec11copy->norm2();
    expectedNorm += vec12copy->Norm2() * vec12copy->Norm2();
    expectedNorm = std::sqrt(expectedNorm);

    bvec1->add(bvec2);
    std::cout << expectedNorm << std::endl << std::flush;
    std::cout << bvec1->norm2() << std::endl << std::flush;
    if (std::abs(expectedNorm-bvec1->norm2()) > TZERO)
        return FAILURE;

    // we also check if adding a non-block vector raises an exception
    try
    {
        bvec1->add(wrap(vec11));
    }
    catch (Exception* e)
    {
        return SUCCESS;
    }

    return FAILURE;
}

// non consistent dimensions
int add2()
{
    unsigned int N = 3;

    shp<BlockVector> bvec1(new BlockVector(2));
    shp<VECTOREPETRA> vec11 = randvecepetra(N);
    shp<DENSEVECTOR> vec12 = randvecdense(N);
    bvec1->setBlock(0, wrap(vec11));
    bvec1->setBlock(1, wrap(vec12));

    shp<VECTOREPETRA> vec11copy(new VECTOREPETRA(*vec11));
    shp<DENSEVECTOR> vec12copy(new DENSEVECTOR(*vec12));

    shp<BlockVector> bvec2(new BlockVector(1));
    shp<VECTOREPETRA> vec21 = randvecepetra(N);
    bvec2->setBlock(0, wrap(vec21));

    try
    {
        bvec1->add(bvec2);
    }
    catch (Exception* e)
    {
        return SUCCESS;
    }

    return FAILURE;
}

// current vector is empty
int add3()
{
    bool areEqual;

    unsigned int M = 4;

    shp<BlockVector> bvec1(new BlockVector(0));
    shp<BlockVector> bvec2(new BlockVector(2));
    bvec2->setBlock(0, wrap(randvecepetra(M)));
    bvec2->setBlock(1, wrap(randvecdense(M+1)));

    bvec1->add(bvec2);

    if (std::abs(bvec1->norm2() - bvec2->norm2()) < TZERO)
        return SUCCESS;

    return FAILURE;
}

// other vector is empty
int add4()
{
    bool areEqual;

    unsigned int M = 5;

    shp<BlockVector> bvec1(new BlockVector(2));
    shp<BlockVector> bvec2(new BlockVector(0));
    bvec1->setBlock(0, wrap(randvecepetra(M)));
    bvec1->setBlock(1, wrap(randvecdense(M+1)));

    double prevnnorm = bvec1->norm2();

    bvec1->add(bvec2);

    if (std::abs(bvec1->norm2() - prevnnorm) < TZERO)
        return SUCCESS;

    return FAILURE;
}

int multiplyByScalar()
{
    unsigned int M = 5;

    shp<BlockVector> bvec1(new BlockVector(2));
    bvec1->setBlock(0, wrap(randvecepetra(M)));
    bvec1->setBlock(1, wrap(randvecdense(M+1)));

    double normm = bvec1->norm2();
    bvec1->multiplyByScalar(1.234);

    if (abs(normm * 1.234 - bvec1->norm2()) < TZERO)
        return SUCCESS;

    return FAILURE;
}

// we just check that no errors are thrown
int _clone()
{
    unsigned int M = 5;

    shp<BlockVector> bvec(new BlockVector(2));
    bvec->setBlock(0, wrap(randvecepetra(M)));
    bvec->setBlock(1, wrap(randvecdense(M+1)));

    shp<aVector> aux(bvec->clone());

    return SUCCESS;
}

int shallowCopy()
{
    unsigned int M = 6;

    shp<BlockVector> bvec1(new BlockVector(2));
    bvec1->setBlock(0, wrap(randvecepetra(M)));
    bvec1->setBlock(1, wrap(randvecdense(M+1)));

    shp<BlockVector> bvec2(new BlockVector(0));
    bvec2->shallowCopy(bvec1);

    bvec2->multiplyByScalar(2);

    // check if also dmat was modified
    if (std::abs(bvec2->norm2() - bvec1->norm2()) < TZERO)
        return SUCCESS;

    return FAILURE;
}

int deepCopy()
{
    unsigned int M = 6;

    shp<BlockVector> bvec1(new BlockVector(2));
    bvec1->setBlock(0, wrap(randvecepetra(M)));
    bvec1->setBlock(1, wrap(randvecdense(M+1)));

    shp<BlockVector> bvec2(new BlockVector(0));
    bvec2->deepCopy(bvec1);

    bvec2->multiplyByScalar(2);

    // check if also dmat was modified
    if (std::abs(bvec2->norm2() - bvec1->norm2()) > TZERO)
        return SUCCESS;

    return FAILURE;
}

int globalTypeIs()
{
    unsigned int M = 6;

    shp<BlockVector> bvec1(new BlockVector(2));
    bvec1->setBlock(0, wrap(randvecepetra(M)));
    bvec1->setBlock(1, wrap(randvecdense(M+1)));

    if (bvec1->globalTypeIs(DISTRIBUTED))
        return FAILURE;

    bvec1->setBlock(1, wrap(randvecepetra(M+1)));
    if (bvec1->globalTypeIs(DISTRIBUTED))
        return SUCCESS;

    return FAILURE;
}

int getSubvector()
{
    unsigned int M = 6;
    shp<BlockVector> bvec(new BlockVector(6));
    bvec->setBlock(0, wrap(randvecepetra(M)));
    bvec->setBlock(1, wrap(randvecdense(M+1)));
    bvec->setBlock(2, wrap(randvecdense(M-1)));
    auto v1 = randvecepetra(M);
    bvec->setBlock(3, wrap(v1));
    auto v2 = randvecdense(M+2);
    bvec->setBlock(4, wrap(v2));
    bvec->setBlock(5, wrap(randvecdense(M-1)));

    auto subvec = bvec->getSubvector(3,4);

    if (subvec->nRows() != 2)
        return FAILURE;

    if (abs(subvec->block(0)->norm2() - v1->norm2()) > TZERO)
        return FAILURE;

    return SUCCESS;
}

int block()
{
    shp<BlockVector> bvec(new BlockVector(10));

    try
    {
        bvec->block(100);
    }
    catch (Exception* e)
    {
        return SUCCESS;
    }
    return FAILURE;
}

int main()
{
    MPI_Init(nullptr, nullptr);
    srand(time(NULL));

    int status = SUCCESS;

    status |= AtomicTest("Constructor1", &constructor1).run();
    status |= AtomicTest("Constructor2", &constructor2).run();
    status |= AtomicTest("Add1", &add1).run();
    status |= AtomicTest("Add2", &add2).run();
    status |= AtomicTest("Add3", &add3).run();
    status |= AtomicTest("Add4", &add4).run();
    status |= AtomicTest("MultiplyByScalar", &multiplyByScalar).run();
    status |= AtomicTest("clone", &_clone).run();
    status |= AtomicTest("ShallowCopy", &shallowCopy).run();
    status |= AtomicTest("DeepCopy", &deepCopy).run();
    status |= AtomicTest("GlobalTypeIs", &globalTypeIs).run();
    status |= AtomicTest("GetSubvector", &getSubvector).run();
    status |= AtomicTest("Block", &block).run();

    return status;
}
