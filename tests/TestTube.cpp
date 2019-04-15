#ifdef EPETRA_MPI
#include <mpi.h>
#include <Epetra_MpiComm.h>
#else
#include <Epetra_SerialComm.h>
#endif

#include <iostream>
#include <string>
#include <Test.hpp>
#include <Tube.hpp>

using namespace ReMA;

void subTest1(Test& test)
{
    Tube tube(test.getComm());
    try
    {
        tube.setParameterValue("abab", 2.0);
        test.assertTrue(false);
    }
    catch (std::exception& e)
    {
        test.assertTrue(true);
    }
}

void subTest2(Test& test)
{
    Tube tube(test.getComm());
    test.assertTrue(tube.readMesh("../geometries/") == 0);
}

int main()
{
    #ifdef HAVE_MPI
    MPI_Init (nullptr, nullptr);
    std::shared_ptr<Epetra_Comm> comm(new Epetra_MpiComm(MPI_COMM_WORLD));
    #else
    std::shared_ptr<Epetra_Comm> comm(new Epetra_SerialComm ());
    #endif

    Test test("TubeTest",comm);
    test.addSubTest(*subTest1);
    test.addSubTest(*subTest2);
    test.run();

    return 0;
}
