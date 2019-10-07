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
#include <mpi.h>
#include <Epetra_MpiComm.h>
#else
#include <Epetra_SerialComm.h>
#endif

#include <functional>

#include <GlobalSolver.hpp>
#include <NavierStokesAssembler.hpp>
#include <PseudoFSIAssembler.hpp>

#include <lifev/core/filter/GetPot.hpp>

#include "ud_functions.hpp"

using namespace RedMA;

double maxLaw(double t)
{
    return LifeV::pressureWomerTime(t,0,0,0,0);
}

double maxLawDt(double t)
{
    return LifeV::pressuredtWomerTime(t,0,0,0,0);
}

int main(int argc, char **argv)
{
    #ifdef HAVE_MPI
    MPI_Init (nullptr, nullptr);
    std::shared_ptr<Epetra_Comm> comm (new Epetra_MpiComm(MPI_COMM_WORLD));
    #else
    std::shared_ptr<Epetra_Comm> comm(new Epetra_SerialComm());
    #endif

    GetPot datafile("data");
    LifeV::Womersley::setParamsFromGetPot(datafile);

    bool verbose = comm->MyPID() == 0;
    GlobalSolver<NavierStokesAssembler> gs(datafile, comm, verbose);
    gs.setExportNorms("norms_nonconforming.txt");
    gs.setLawInflow(std::function<double(double)>(maxLaw));
    gs.setLawDtInflow(std::function<double(double)>(maxLawDt));
    gs.solve();

    return 0;
}
