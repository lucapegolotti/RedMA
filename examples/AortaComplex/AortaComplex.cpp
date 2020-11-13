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

#include <Epetra_ConfigDefs.h>
#ifdef EPETRA_MPI
#include <mpi.h>
#include <Epetra_MpiComm.h>
#else
#include <Epetra_SerialComm.h>
#endif

#include <redma/RedMA.hpp>
#include <redma/problem/ProblemFEM.hpp>
#include <redma/problem/DataContainer.hpp>

using namespace RedMA;

double distalPressure0(double t)
{
    return 0;
}

double distalPressure1(double t)
{
    return 0;
}

double distalPressure2(double t)
{
    return 0;
}

double distalPressure3(double t)
{
    return 0;
}

double distalPressure4(double t)
{
    return 0;
}

int main(int argc, char **argv)
{
    #ifdef HAVE_MPI
    MPI_Init (nullptr, nullptr);
    std::shared_ptr<Epetra_Comm> comm (new Epetra_MpiComm(MPI_COMM_WORLD));
    #else
    std::shared_ptr<Epetra_Comm> comm(new Epetra_SerialComm ());
    #endif

    DataContainer data;
    data.setDatafile("datafiles/data");
    data.setVerbose(comm->MyPID() == 0);
    data.setDistalPressure(distalPressure0, 0);
    data.setDistalPressure(distalPressure1, 1);
    data.setDistalPressure(distalPressure2, 2);
    data.setDistalPressure(distalPressure3, 3);
    data.setDistalPressure(distalPressure4, 4);
    data.finalize();

    ProblemFEM femProblem(data, comm);

    femProblem.solve();

    return 0;
}
