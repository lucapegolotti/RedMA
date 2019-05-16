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
#include <Test.hpp>
#include <BifurcationSymmetric.hpp>

using namespace RedMA;

void subTest1(Test& test)
{
    BifurcationSymmetric bifurcation(test.getComm());
    try
    {
        bifurcation.setParameterValue("abab", 2.0);
        test.assertTrue(false);
    }
    catch (std::exception& e)
    {
        test.assertTrue(true);
    }
}

void subTest2(Test& test)
{
    BifurcationSymmetric bifurcation(test.getComm());
    test.assertTrue(bifurcation.readMesh("../../meshes/") == 0);
}

int main()
{
    #ifdef HAVE_MPI
    MPI_Init (nullptr, nullptr);
    std::shared_ptr<Epetra_Comm> comm(new Epetra_MpiComm(MPI_COMM_WORLD));
    #else
    std::shared_ptr<Epetra_Comm> comm(new Epetra_SerialComm ());
    #endif

    Test test("BifurcationSymmetric",comm);
    test.addSubTest(*subTest1);
    test.addSubTest(*subTest2);
    test.run();

    return 0;
}
