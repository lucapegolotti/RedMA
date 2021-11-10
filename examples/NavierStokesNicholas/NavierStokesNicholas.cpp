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

#include <redma/reduced_basis/MatricesGenerator.hpp>
#include <redma/reduced_basis/SnapshotsSampler.hpp>

using namespace RedMA;

double inflow(double t, double a, double c)
{
    double T = 0.3;
    return 1-cos(2*M_PI*t/T) + c*sin(2*M_PI*a*t/T);
}

int main(int argc, char **argv)
{
    
    #ifdef HAVE_MPI
    MPI_Init (nullptr, nullptr);
    EPETRACOMM comm (new Epetra_MpiComm(MPI_COMM_WORLD));
    #else
    EPETRACOMM comm(new Epetra_SerialComm());
    #endif
    printlog(MAGENTA,"start example", true);	
    DataContainer data;
    data.setDatafile("datafiles/data_fem");
    data.setVerbose(comm->MyPID() == 0);

    SnapshotsSampler sampler(data, inflow, comm);
    sampler.takeSnapshots();

    return 0;
}
