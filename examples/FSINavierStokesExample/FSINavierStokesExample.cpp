//
// Created by micol on 13/03/2021.
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
    return (t<0.02)*(1.0 - std::cos((t) * M_PI / (0.02 )) )/ 2.0 +(t>=0.02);
    //return (t<0.02)*t+ (t>0.02)*0.02 ;
     //return 1;
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
    data.setDatafile("datafiles/data_fem");
    //data.setInflow(inflow);
    data.setVerbose(comm->MyPID() == 0);
    data.finalize();
    GlobalProblem femProblem(data, comm);

    femProblem.solve();

    return 0;
}
