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

#include <Epetra_ConfigDefs.h>
#ifdef EPETRA_MPI
//#include <mpi.h>
#include <Epetra_MpiComm.h>
#else
#include <Epetra_SerialComm.h>
#endif

#include <lifev/core/LifeV.hpp>
#include <lifev/core/util/LifeChronoManager.hpp>

#include <rb/reduced_basis/rbLifeV.hpp>

#include <lifev/navier_stokes_blocks/solver/NavierStokesOperator.hpp>
#include <lifev/core/linear_algebra/ApproximatedInvertibleRowMatrix.hpp>
#include <lifev/navier_stokes_blocks/solver/NavierStokesPreconditionerOperator.hpp>
#include <lifev/navier_stokes_blocks/solver/aSIMPLEOperator.hpp>

// using namespace LifeV;
// using namespace rbLifeV;

int main (int argc, char** argv)
{
    // MPI initialization
#ifdef HAVE_MPI
    MPI_Init(&argc, &argv);
    std::shared_ptr<Epetra_Comm> Comm(new Epetra_MpiComm(MPI_COMM_WORLD));
#else
    std::shared_ptr<Epetra_Comm> Comm(new Epetra_SerialComm);
#endif


    return 0;
}
