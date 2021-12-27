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
#include <chrono>
#include <thread>
#include <time.h>       /* time */

using namespace RedMA;

int main(int argc, char **argv)
{
    srand(1234);

    #ifdef HAVE_MPI
    MPI_Init (nullptr, nullptr);
    EPETRACOMM comm (new Epetra_MpiComm(MPI_COMM_WORLD));
    #else
    EPETRACOMM comm(new Epetra_SerialComm());
    #endif

    Chrono chrono;
    chrono.start();

    std::string msg = "Starting chrono\n";
    printlog(MAGENTA, msg, true);

    DataContainer data;
    data.setDatafile("datafiles/data_rb");
    data.setVerbose(comm->MyPID() == 0);

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

    GlobalProblem rbProblem(data, comm);
    rbProblem.solve();

    msg = "Saved the solution at " + curSolutionDir + "\n";
    printlog(MAGENTA, msg, true);

    msg = "Total time =  ";
    msg += std::to_string(chrono.diff());
    msg += " seconds\n";
    printlog(MAGENTA, msg, true);

    return 0;
}
