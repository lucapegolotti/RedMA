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

using namespace RedMA;

double inletDirichlet0(double t)
{
    return 1;
}

double inletDirichlet1(double t)
{
    return 1;
}

double inletNeumann0(double t)
{
    const double t0 = 0.0;
    const double T = t0 + 3e-3;
    const double omega = 2.0 * M_PI / T;
    const double Pmax = 13300.0;

    if ((t >= t0) && (t <= T))
        return -0.5 * (1.0 - std::cos(omega * (t-t0))) * Pmax;

    return 0.0;
}

double inletNeumann1(double t)
{
    const double t0 = 1e-2;
    const double T = t0 + 3e-3;
    const double omega = 2.0 * M_PI / T;
    const double Pmax = 13300.0;

    if ((t >= t0) && (t <= T))
        return -0.5 * (1.0 - std::cos(omega * (t-t0))) * Pmax;

    return 0.0;
}

int main(int argc, char **argv)
{
    #ifdef HAVE_MPI
    MPI_Init (nullptr, nullptr);
    EPETRACOMM comm (new Epetra_MpiComm(MPI_COMM_WORLD));
    #else
    EPETRACOMM comm(new Epetra_SerialComm());
    #endif

    Chrono chrono;
    chrono.start();

    std::string msg = "Starting chrono \n";
    printlog(MAGENTA, msg, true);

    DataContainer data;
    data.setDatafile("datafiles/data");
    data.setVerbose(comm->MyPID() == 0);

    if (!std::strcmp(data("bc_conditions/inlet_bc_type", "dirichlet").c_str(), "dirichlet"))
    {
        // the second argument is teh index of the inlet. The corresponding flag is set in the datafile
        data.setInletBC(inletDirichlet0, 0);
        data.setInletBC(inletDirichlet1, 1);
    }
    else if (!std::strcmp(data("bc_conditions/inlet_bc_type", "dirichlet").c_str(), "neumann"))
    {
        // the second argument is teh index of the inlet. The corresponding flag is set in the datafile
        data.setInletBC(inletNeumann0, 0);
        data.setInletBC(inletNeumann1, 1);
    }
    else
        throw new Exception("Unrecognized inlet BC type!");

    std::string solutionDir = "solutions";
    if (!fs::exists(solutionDir))
        fs::create_directory(solutionDir);

    unsigned int solutionIndex = 0;
    std::string curSolutionDir = solutionDir + "/solution" + std::to_string(solutionIndex) + "/";
    while (fs::exists(curSolutionDir))
    {
        solutionIndex++;
        curSolutionDir = solutionDir + "/solution" + std::to_string(solutionIndex) + "/";
    }

    if (solutionIndex >= 100){
        msg = "Cannot store more than 100 solutions!\n";
        printlog(RED, msg, true);
        exit(1);
    }
    else {
        fs::create_directory(curSolutionDir);
        msg = "Saving the solution at " + curSolutionDir + "\n";
        printlog(MAGENTA, msg, true);
    }

    data.setValueString("exporter/outdir", curSolutionDir);
    data.finalize();

    GlobalProblem globalProblem(data, comm);
    globalProblem.solve();

    msg = "Saved the solution at " + curSolutionDir + "\n";
    printlog(MAGENTA, msg, true);

    msg = "Total time =  ";
    msg += std::to_string(chrono.diff());
    msg += " seconds\n";
    printlog(MAGENTA, msg, true);

    GlobalProblem femProblem(data, comm);

    femProblem.solve();

    return 0;
}
