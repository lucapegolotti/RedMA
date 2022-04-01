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
#include <redma/problem/DataContainer.hpp>
#include <redma/problem/GlobalProblem.hpp>
#include <thread>

using namespace RedMA;

double inflow(double t)
{
    double T = 0.3;
    return 1-cos(2*M_PI*t/T) + 0.2*sin(2*6*M_PI*t/T);
}

int main(int argc, char **argv)
{
    srand(1234);

    #ifdef HAVE_MPI
    MPI_Init (nullptr, nullptr);
    EPETRACOMM comm (new Epetra_MpiComm(MPI_COMM_WORLD));
    #else
    EPETRACOMM comm(new Epetra_SerialComm());
    #endif

    std::string msg = "Starting RB matrices generation\n";
    printlog(MAGENTA, msg, true);

    DataContainer data;
    data.setDatafile("datafiles/data_rb");
    data.setVerbose(comm->MyPID() == 0);
    data.setInletBC(inflow);
    data.finalize();

    GlobalProblem rbProblem(data, comm);
    // rbProblem.solve();

    msg = "Done\n";
    printlog(MAGENTA, msg, true);

    return 0;
}