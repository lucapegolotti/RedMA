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

#include <thread>

#include "inflows.hpp"

#include <redma/reduced_basis/SnapshotsSampler.hpp>

using namespace RedMA;

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

    double T = data("time_discretization/T", 1.0);
    double Tramp = - data("time_discretization/t0ramp", 0.05);

    /*if (std::strcmp(data("rb/offline/snapshots/param_type", "inflow").c_str(), "inflow"))
        throw new Exception("This test case handles only 'inflow' parametrization!");*/

    SnapshotsSampler sampler(data, comm);
    if (!std::strcmp(data("rb/offline/snapshots/inflow_type", "default").c_str(), "default"))
        sampler.setInflow([T](const double t, const std::vector<double> params){return inflow(t, params, T);});
    else if (!std::strcmp(data("rb/offline/snapshots/inflow_type", "default").c_str(), "periodic"))
        sampler.setInflow([T, Tramp](const double t, const std::vector<double> params){return inflow_periodic(t, params, T, Tramp);});
    else if (!std::strcmp(data("rb/offline/snapshots/inflow_type", "default").c_str(), "systolic"))
        sampler.setInflow([Tramp](const double t, const std::vector<double> params){return inflow_systolic(t, params, Tramp);});
    else if (!std::strcmp(data("rb/offline/snapshots/inflow_type", "default").c_str(), "heartbeat"))
        sampler.setInflow([Tramp](const double t, const std::vector<double> params){return inflow_heartbeat(t, params, Tramp);});
    else if (!std::strcmp(data("rb/offline/snapshots/inflow_type", "default").c_str(), "bypass"))
        sampler.setInflow([Tramp](const double t, const std::vector<double> params){return inflow_bypass(t, params, Tramp);});
    else
        throw new Exception("Unrecognized type of inflow parametrization! "
                            "Available types: {default, periodic, systolic, heartbeat}.");

    sampler.takeSnapshots(Nstart);

    // To solve the RB problem in RedMA  --> set datafile accordingly!
    /*data.setInletBC(inflow2, 0);
    data.finalize();
    GlobalProblem rbProblem(data, comm);
    rbProblem.solve();*/

    // To compute the convective term, given a solution
    /*GlobalProblem femProblem(data, comm);
    femProblem.validateRBConvectiveTerm();*/

    msg = "Total time =  ";
    msg += std::to_string(chrono.diff());
    msg += " seconds\n";
    printlog(MAGENTA, msg, true);

    return 0;
}
