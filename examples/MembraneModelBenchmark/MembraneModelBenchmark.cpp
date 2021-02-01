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

using namespace RedMA;

/*#define COEFF 10*/

/*double maxLaw(double t)
{
    const double T = 3e-3;
    const double omega = 2.0 * M_PI / (T);
    const double Pmax = 13300.0;
    if (t <= T)
    {
        return -0.5 * (1.0 - cos(omega * t) ) * Pmax;
    }
    return 0;
}

double maxLawDt(double t)
{
    const double T = 3e-3;
    const double omega = 2.0 * M_PI / (T);
    const double Pmax = 13300.0;
    if (t <= T)
    {
        return -0.5 * omega * sin(omega * t) * Pmax;
    }
    return 0.0;
}*/

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

    std::string msg = "Starting chrono...\n";
    printlog(MAGENTA, msg, true);

    DataContainer data;
    data.setDatafile("datafiles/data_fem");
    data.setVerbose(comm->MyPID() == 0);
    std::string outdir = "solution_fem_reference/";
    data.setValueString("exporter/outdir", outdir);
    data.finalize();

    /*gs.setExportNorms("norms_nonconforming.txt");
    gs.setLawInflow(std::function<double(double)>(maxLaw));
    gs.setLawDtInflow(std::function<double(double)>(maxLawDt));*/

    GlobalProblem globalProblem(data, comm);
    globalProblem.solve();

    msg = "Total time =  ";
    msg += std::to_string(chrono.diff());
    msg += " seconds\n";
    printlog(MAGENTA, msg, true);

    return 0;
}
