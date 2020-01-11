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

#include <stdlib.h>
#include <time.h>
#include <iostream>
#include <string>

#include <redma/utils/Test.hpp>
#include <redma/geometry/Tube.hpp>

using namespace RedMA;

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
    test.assertTrue(tube.readMesh("../../meshes/") == 0);
}

void subTest3(Test& test)
{
    typedef std::shared_ptr<Tube> TubePtr;
    typedef LifeV::VectorSmall<3> Vector3D;
    // we impose a random affine transformation and check if the child tube
    // stay attached to the first one
    srand(time(NULL));

    TubePtr tube1(new Tube(test.getComm()));
    TubePtr tube2(new Tube(test.getComm()));
    tube2->setIsChild(true);

    tube1->setIsChild(false);
    tube1->setParameterValue("alphax", static_cast<float>(rand()) / RAND_MAX);
    tube1->setParameterValue("alphay", static_cast<float>(rand()) / RAND_MAX);
    tube1->setParameterValue("alphaz", static_cast<float>(rand()) / RAND_MAX);
    tube1->setParameterValue("bx", rand());
    tube1->setParameterValue("by", rand());
    tube1->setParameterValue("bz", rand());
    tube1->setParameterValue("scale", 0.888);
    try
    {
        // we check if assert on mesh being read holds
        tube1->applyAffineTransformation();
        test.assertTrue(false);
    }
    catch (Exception& e)
    {
        test.assertTrue(true);
        tube1->readMesh();
        tube1->applyAffineTransformation();
    }
    GeometricFace inlet = tube1->getInlet();

    const double tol = 1e-14;

    // we check if the inlet has been moved
    Vector3D res = inlet.M_center - Vector3D(0,0,0);
    test.assertTrue(res.norm() > tol);

    res = inlet.M_normal - Vector3D(0,0,-1);
    test.assertTrue(res.norm() > tol);

    tube2->mapChildInletToParentOutlet(tube1->getOutlet(0));
    tube2->readMesh();
    tube2->applyAffineTransformation();

    res = tube1->getOutlet(0).M_center - tube2->getInlet().M_center;
    test.assertTrue(res.norm() < tol);

    // normals must have opposite sign
    res = tube1->getOutlet(0).M_normal + tube2->getInlet().M_normal;
    test.assertTrue(res.norm() < tol);

    test.assertTrue(std::abs(tube1->getOutlet(0).M_radius -
                             tube2->getInlet().M_radius) < tol);
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
    test.addSubTest(*subTest3);
    test.run();

    return 0;
}
