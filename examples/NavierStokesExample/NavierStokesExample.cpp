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

using namespace RedMA;

#define COEFF 4
#define HEARTBEAT 0.8
// double maxLaw(double t)
// {
//     double tt = std::fmod(t,HEARTBEAT);
//     double poly = 0;
//     const double coeffs[8] = {-743.0, 2524.6, -3318.4, 2087.0, -606.25, 49.5, 6.6, 0.0};
//
//     double mon = 1.0;
//     for (int i = 0; i < 8; i++)
//     {
//         poly += coeffs[7-i] * mon;
//         mon *= tt / HEARTBEAT;
//     }
//     return poly * COEFF;
// }
//
// double maxLawDt(double t)
// {
//     double tt = std::fmod(t,HEARTBEAT);
//     double poly = 0;
//     const double coeffs[7] = {-743.0, 2524.6, -3318.4, 2087.0, -606.25, 49.5, 6.6};
//
//     double mon = 1.0;
//     for (int i = 0; i < 7; i++)
//     {
//         poly += coeffs[6-i] * mon;
//         mon *= tt / HEARTBEAT;
//     }
//     return poly * COEFF;
// }

double maxLaw_(double t)
{
    return std::sin(t);
}

double maxLawDt_(double t)
{
    return std::cos(t);
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

    bool verbose = comm->MyPID() == 0;

    GlobalSolver<NavierStokesAssembler> gs(datafile, comm, verbose);
    gs.setExportNorms("norms_nonconforming.txt");
    gs.setLawInflow(std::function<double(double)>(maxLaw_));
    gs.setLawDtInflow(std::function<double(double)>(maxLawDt_));
    gs.solve();

    return 0;
}
