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

double inletDirichlet(double t)
{
    const double T = 5e-3;
    const double omega = M_PI / T;
    const double Q_max = 1.0;

    if (t <= T)
        return 0.5 * (1.0 - std::cos(omega * t)) * Q_max;

    return Q_max;
}

double inletNeumann(double t)
{
    const double T = 3e-3;
    const double omega = 2.0 * M_PI / T;
    const double P_max = 13300.0;

    if (t <= T)
        return -0.5 * (1.0 - std::cos(omega * t) ) * P_max;

    return 0.0;
}

double outletNeumann(double t, unsigned int i)
{
    const double T = 5e-3;
    const double omega = M_PI / T;
    double P_max;
    // potentially we can define here different BC for different outlet indices
    if (i == 0)
        P_max = 1.5 * 1333.0;
    else
        P_max = 1333.0;

    if (t <= T)
        return 0.5 * (1.0 - std::cos(omega * t)) * P_max;

    return P_max;
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

    std::string msg = "Starting chrono... \n";
    printlog(MAGENTA, msg, true);

    DataContainer data;
    data.setDatafile("datafiles/data");
    data.setVerbose(comm->MyPID() == 0);

    if (!std::strcmp(data("bc_conditions/inlet_bc_type", "dirichlet").c_str(), "dirichlet"))
        data.setInletBC(inletDirichlet);
    else if (!std::strcmp(data("bc_conditions/inlet_bc_type", "dirichlet").c_str(), "neumann"))
        data.setInletBC(inletNeumann);
    else
        throw new Exception("Unrecognized inlet BC type! "
                            "Available types: {dirichlet, neumann}.");

    unsigned int numOutletConditions = data("bc_conditions/numoutletbcs", 0);
    for (unsigned int i = 0; i < numOutletConditions; i++)
    {
        std::string dataEntry = "bc_conditions/outlet" + std::to_string(i);
        if (!std::strcmp(data(dataEntry + "/type", "windkessel").c_str(), "neumann"))
            data.setOutletBC([i](double t){return outletNeumann(t,i);}, i);
    }

    std::string solutionDir = "solutions";
    if (!fs::exists(solutionDir))
        fs::create_directory(solutionDir);

    unsigned int solutionIndex = 0;
    std::string curSolutionDir = solutionDir + "/solution" + std::to_string(solutionIndex) + "/";
    while ((fs::exists(curSolutionDir)) && (solutionIndex < 1000))
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

    GlobalProblem femProblem(data, comm);
    femProblem.solve();

    msg = "Total time =  ";
    msg += std::to_string(chrono.diff());
    msg += " seconds\n";
    printlog(MAGENTA, msg, true);

    return 0;
}
