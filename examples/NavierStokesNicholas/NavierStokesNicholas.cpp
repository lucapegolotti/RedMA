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
#include <redma/problem/GlobalProblem.hpp>
#include <redma/problem/DataContainer.hpp>

#include <cmath>
#include <thread>

#include <redma/reduced_basis/SnapshotsSampler.hpp>

using namespace RedMA;

double inflow(double t, double a, double c)
{
    double T = 0.3;
    return 1-cos(2*M_PI*t/T) + c*sin(2*M_PI*a*t/T);
}

int main(int argc, char **argv)
{

    std::mt19937_64 eng{std::random_device{}()};
    std::uniform_int_distribution<> dist{1, 20};
    std::this_thread::sleep_for(std::chrono::seconds{dist(eng)});

    Chrono chrono;
    chrono.start();

    std::string msg = "Starting chrono... \n";
    printlog(MAGENTA, msg, true);
    
    #ifdef HAVE_MPI
    MPI_Init (nullptr, nullptr);
    EPETRACOMM comm (new Epetra_MpiComm(MPI_COMM_WORLD));
    #else
    EPETRACOMM comm(new Epetra_SerialComm());
    #endif

    printlog(MAGENTA,"Starting snapshots generation", true);
    DataContainer data;
    data.setDatafile("datafiles/data_fem");
    data.setVerbose(comm->MyPID() == 0);

    unsigned int Nstart = 0;
    if (argc > 1)
        Nstart = std::atoi(argv[1]);

    SnapshotsSampler sampler(data, comm);
    sampler.setInflow(inflow);
    sampler.takeSnapshots(Nstart);

    msg = "Total time =  ";
    msg += std::to_string(chrono.diff());
    msg += " seconds\n";
    printlog(MAGENTA, msg, true);

    return 0;
}
