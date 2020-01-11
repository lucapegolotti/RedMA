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

#ifdef EPETRA_MPI
#include <mpi.h>
#include <Epetra_MpiComm.h>
#else
#include <Epetra_SerialComm.h>
#endif

#include <iostream>
#include <string>
#include <redma/utils/Test.hpp>
#include <redma/geometry/GeometryParser.hpp>

using namespace RedMA;

void subTest1(Test& test)
{
    try
    {
        // test if trying to open non existing file raises exception
        GeometryParser gParser("Idontexist.xml", test.getComm(), false);
        test.assertTrue(false);
    }
    catch (Exception& e)
    {
        test.assertTrue(true);
    }
}

void subTest2(Test& test)
{
    GeometryParser gParser("testdata/testdata1.xml", test.getComm(), false);

    // we test if the tree has been filled
    test.assertTrue(!gParser.getTree().isEmpty());
}

int main()
{
    #ifdef HAVE_MPI
    MPI_Init (nullptr, nullptr);
    std::shared_ptr<Epetra_Comm> comm(new Epetra_MpiComm(MPI_COMM_WORLD));
    #else
    std::shared_ptr<Epetra_Comm> comm(new Epetra_SerialComm ());
    #endif

    Test test("TubeGeometryParser",comm);
    test.addSubTest(*subTest1);
    test.addSubTest(*subTest2);
    test.run();

    return 0;
}
