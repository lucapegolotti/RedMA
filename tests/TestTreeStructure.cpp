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
#include <BifurcationSymmetric.hpp>
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
    try
    {
        tree.addChild(0,tubePtr2);
        test.assertTrue(true);
    }
    catch(Exception& e)
    {
        test.assertTrue(false);
    }

    std::shared_ptr<Tube> tubePtr3(new Tube(test.getComm()));
    try
    {
        // this fails because tube can only have one child
        tree.addChild(0,tubePtr3);
        test.assertTrue(false);
    }
    catch(Exception& e)
    {
        test.assertTrue(true);
    }
}

void subTest4(Test& test)
{
    TreeStructure tree;
    std::shared_ptr<Tube> tubePtr1(new Tube(test.getComm()));
    tree.setRoot(tubePtr1);

    std::shared_ptr<BifurcationSymmetric>
                            bifurcationPtr2(new BifurcationSymmetric(test.getComm()));

    unsigned int id = id = tree.addChild(0,bifurcationPtr2);

    std::shared_ptr<Tube> tubePtr3(new Tube(test.getComm()));
    tree.addChild(id,tubePtr3);

    std::shared_ptr<Tube> tubePtr5(new Tube(test.getComm()));
    // bifurcation should allow for two children
    try
    {
        tree.addChild(id,tubePtr3);
        test.assertTrue(true);
    }
    catch(Exception& e)
    {
        test.assertTrue(false);
    }

    test.assertTrue(tree.depth() == 2);
}

int main()
{
    #ifdef HAVE_MPI
    MPI_Init (nullptr, nullptr);
    std::shared_ptr<Epetra_Comm> comm(new Epetra_MpiComm(MPI_COMM_WORLD));
    #else
    std::shared_ptr<Epetra_Comm> comm(new Epetra_SerialComm());
    #endif

    Test test("TreeStructureTest",comm);
    test.addSubTest(*subTest1);
    test.addSubTest(*subTest2);
    test.addSubTest(*subTest3);
    test.addSubTest(*subTest4);
    test.run();

    return 0;
}
