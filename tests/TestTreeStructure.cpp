// Reduced Modeling of Arteries (ReMA)
// Copyright (C) 2019  Luca Pegolotti
//
// ReMA is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// ReMA is distributed in the hope that it will be useful,
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

#include <TreeStructure.hpp>
#include <Tube.hpp>
#include <Test.hpp>
#include <Exception.hpp>

using namespace ReMA;

void subTest1(Test& test)
{
    TreeStructure tree;
    test.assertTrue(tree.getMaxID() == 0);
}

void subTest2(Test& test)
{
    TreeStructure tree;
    std::shared_ptr<Tube> tubePtr(new Tube(test.getComm()));
    try
    {
        tree.setRoot(tubePtr);
        test.assertTrue(true);
    }
    catch(Exception& e)
    {
        test.assertTrue(false);
    }
}

void subTest3(Test& test)
{
    TreeStructure tree;
    std::shared_ptr<Tube> tubePtr1(new Tube(test.getComm()));
    tree.setRoot(tubePtr1);

    std::shared_ptr<Tube> tubePtr2(new Tube(test.getComm()));
    tree.addChild(0,tubePtr2);
}

int main()
{
    #ifdef HAVE_MPI
    MPI_Init (nullptr, nullptr);
    std::shared_ptr<Epetra_Comm> comm(new Epetra_MpiComm(MPI_COMM_WORLD));
    #else
    std::shared_ptr<Epetra_Comm> comm(new Epetra_SerialComm ());
    #endif

    Test test("TreeStructureTest",comm);
    test.addSubTest(*subTest1);
    test.addSubTest(*subTest2);
    test.addSubTest(*subTest3);
    test.run();

    return 0;
}
