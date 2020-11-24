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
#include <redma/array/DenseVector.hpp>

#include <tests/AtomicTest.hpp>

#include <time.h>

using namespace RedMA;

shp<DENSEVECTOR> randvec(unsigned int N)
{
    shp<DENSEVECTOR> vec(new DENSEVECTOR(N));

    for (unsigned int i = 0; i < N; i++)
        vec->operator()(i) = rand() % 10;

    return vec;
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
    DenseVector dVector;

    return SUCCESS;
}

int setVector()
{
    shp<DENSEVECTOR> vec = randvec(5);
    DenseVector dVector;
    dVector.setVector(vec);

    bool areEqual = checkEqual(dVector.getVector(), vec);

    return !areEqual;
}

// case in which both vectors are filled
int add1()
{
    unsigned int N = 3;

    shp<DENSEVECTOR> vec1 = randvec(N);
    shp<DENSEVECTOR> vec1Copy(new DENSEVECTOR(*vec1));
    shp<DENSEVECTOR> vec2 = randvec(N);

    shp<DenseVector> dvec1(new DenseVector());
    dvec1->setVector(vec1);

    shp<DenseVector> dvec2(new DenseVector());
    dvec2->setVector(vec2);

    dvec1->add(dvec2);

    for (unsigned int i = 0; i < N; i++)
    {
        double expected = vec1Copy->operator()(i) + vec2->operator()(i);
        if (abs(dvec1->getVector()->operator()(i) - expected) > TZERO)
            return FAILURE;
    }
    return SUCCESS;
}

// non consistent dimensions
int add2()
{
    unsigned int N = 3;

    shp<DENSEVECTOR> vec1 = randvec(N);
    shp<DENSEVECTOR> vec2 = randvec(N+1);

    shp<DenseVector> dvec1(new DenseVector());
    dvec1->setVector(vec1);

    shp<DenseVector> dvec2(new DenseVector());
    dvec2->setVector(vec2);

    try
    {
        dvec1->add(dvec2);
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

    shp<DenseVector> dvec1(new DenseVector());
    shp<DenseVector> dvec2(new DenseVector());
    shp<DENSEVECTOR> vec2 = randvec(M);

    dvec2->setVector(vec2);
    dvec1->add(dvec2);
    areEqual = checkEqual(vec2, dvec1->getVector());

    if (areEqual)
        return SUCCESS;

    return FAILURE;
}

// other vector is empty
int add4()
{
    bool areEqual;

    unsigned int N = 5;

    shp<DenseVector> dvec1(new DenseVector());
    shp<DenseVector> dvec2(new DenseVector());
    shp<DENSEVECTOR> vec1 = randvec(N);

    dvec1->setVector(vec1);
    dvec1->add(dvec2);

    areEqual = checkEqual(vec1, dvec1->getVector());

    if (areEqual)
        return SUCCESS;

    return FAILURE;
}

int multiplyByScalar()
{
    unsigned int N = 5;

    shp<DenseVector> dvec1(new DenseVector());
    shp<DENSEVECTOR> vec1 = randvec(N);

    double normm = vec1->NormInf();

    dvec1->setVector(vec1);
    dvec1->multiplyByScalar(1.234);

    if (abs(normm * 1.234 - dvec1->getVector()->NormInf()) < TZERO)
        return SUCCESS;

    return FAILURE;
}
// we just check that no errors are thrown
int clone()
{
    unsigned int M = 5;

    shp<DenseVector> dvec(new DenseVector());
    shp<DENSEVECTOR> vec = randvec(M);
    dvec->setVector(vec);

    shp<aVector> aux(dvec->clone());

    return SUCCESS;
}

int shallowCopy()
{
    unsigned int N = 6;

    shp<DenseVector> dvec(new DenseVector());
    dvec->setVector(randvec(N));

    shp<DenseVector> otherdvec(new DenseVector());
    otherdvec->shallowCopy(dvec);

    otherdvec->getVector()->operator()(2) = otherdvec->getVector()->operator()(2) + 1;

    // check if also dmat was modified
    if (otherdvec->getVector()->operator()(2) == dvec->getVector()->operator()(2))
        return SUCCESS;

    return FAILURE;
}

int deepCopy()
{
    unsigned int N = 6;

    shp<DenseVector> dvec(new DenseVector());
    dvec->setVector(randvec(N));

    shp<DenseVector> otherdvec(new DenseVector());
    otherdvec->deepCopy(dvec);

    otherdvec->getVector()->operator()(2) = otherdvec->getVector()->operator()(2) + 1;

    // check if dvec wasn't modified
    if (otherdvec->getVector()->operator()(2) != dvec->getVector()->operator()(2))
        return SUCCESS;

    return FAILURE;
}

int accessOperator()
{
    unsigned int N = 6;

    DenseVector dvec;
    dvec.setVector(randvec(N));

    if (dvec(3) == dvec.getVector()->operator()(3))
        return SUCCESS;

    return FAILURE;
}

int main()
{
    srand(time(NULL));

    int status = SUCCESS;

    status |= AtomicTest("Constructor", &constructor).run();
    status |= AtomicTest("SetVector", &setVector).run();
    status |= AtomicTest("Add1", &add1).run();
    status |= AtomicTest("Add2", &add2).run();
    status |= AtomicTest("Add3", &add3).run();
    status |= AtomicTest("Add4", &add4).run();
    status |= AtomicTest("MultiplyByScalar", &multiplyByScalar).run();
    status |= AtomicTest("clone", &clone).run();
    status |= AtomicTest("ShallowCopy", &shallowCopy).run();
    status |= AtomicTest("DeepCopy", &deepCopy).run();
    status |= AtomicTest("Access", &accessOperator).run();

    return status;
}
