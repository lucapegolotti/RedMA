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
#include <redma/solver/problem/ProblemFEM.hpp>
#include <redma/solver/problem/DataContainer.hpp>
#include <redma/geometry/CommunicatorsDistributor.hpp>

using namespace RedMA;

double inflow(double t)
{
    return 10*sin(t);
}

double inflowDt(double t)
{
    return 10*cos(t);
}

int main(int argc, char **argv)
{
    #ifdef HAVE_MPI
    MPI_Init (nullptr, nullptr);
    EPETRACOMM comm (new Epetra_MpiComm(MPI_COMM_WORLD));
    #else
    EPETRACOMM comm(new Epetra_SerialComm());
    #endif

    bool distributed = true;

    DataContainer data;
    data.setDatafile("datafiles/data");
    data.setInflow(inflow);
    data.setInflowDt(inflowDt);
    data.setMasterComm(comm);
    data.setDistributed(distributed);
    data.setVerbose(comm->MyPID() == 0);

    if (distributed)
    {
        CommunicatorsDistributor commD(data.getDatafile(), comm);
        ProblemFEM femProblem(data, commD);
        femProblem.solve();
    }
    else
    {
        ProblemFEM femProblem(data, comm);
        femProblem.solve();
    }

    return 0;
}
