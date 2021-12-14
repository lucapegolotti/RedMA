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

double inflow(double t)
{
    return t;
}

int main(int argc, char **argv)
{
#ifdef HAVE_MPI
    MPI_Init (nullptr, nullptr);
    EPETRACOMM comm (new Epetra_MpiComm(MPI_COMM_WORLD));
#else
    EPETRACOMM comm(new Epetra_SerialComm());
#endif

    std::string msg = "Starting chrono...\n";
    printlog(MAGENTA, msg, true);

    DataContainer data;
    data.setDatafile("datafiles/data_fem");
    data.setInflow(inflow);
    data.setVerbose(comm->MyPID() == 0);

    GlobalProblem femProblem(data, comm);

    std::string path = data("importer/indir", "");
    int param_nb = data("importer/param_nb", -1);
    bool isFEM = data("importer/isFEM", false);

    if ((param_nb<0) || (!std::strcmp(path.c_str(), "")))
        throw new Exception("Invalid importing path or invalid parameter number!");
    else
        path += ("param" + std::to_string(param_nb) + "/");

    if (isFEM)
         path +=  "FEM/";
    else
        path += "RB/";

    femProblem.exportFromFiles(path);

    return 0;
}


