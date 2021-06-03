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

double inflow2(double t)
{
    //return 1;
    return (t<0.02)*(1.0 - std::cos((t) * M_PI / (0.02 )) )/ 2.0 +(t>=0.02);
}

double inflow3(double t)
{
    //return 1;
    return (t<0.02)*(1.0 - std::cos((t) * M_PI / (0.02 )) )/ 2.0 +(t>=0.02);
}

int main(int argc, char **argv)
{
    #ifdef HAVE_MPI
    MPI_Init (nullptr, nullptr);
    EPETRACOMM comm (new Epetra_MpiComm(MPI_COMM_WORLD));
    #else
    EPETRACOMM comm(new Epetra_SerialComm());
    #endif

    DataContainer data;
    data.setDatafile("datafiles/data_rb");
    data.setVerbose(comm->MyPID() == 0);

    bool analyticInflow = true;
    if (analyticInflow)
    {
        double alpha=0.2125;
        auto inflow1_alpha([alpha](double t) {
            return alpha * inflow2(t);
        });
        auto inflow2_alpha([alpha](double t) {
            return (1 - alpha) * inflow3(t);
        });
        // the second argument is the flag of the inlet. It depends on the mesh.
        data.setInflow(inflow1_alpha, 2);
        data.setInflow(inflow2_alpha, 3);
    }
    else
    {
        data.finalize();
    }

    GlobalProblem femProblem(data, comm);

    femProblem.solve();

    return 0;
}
