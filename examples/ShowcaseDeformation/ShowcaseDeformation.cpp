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

#include <redma/problem/DataContainer.hpp>

// #include <redma/geometry/building_blocks/Tube.hpp>
// #include <redma/geometry/building_blocks/BifurcationSymmetric.hpp>
#include <redma/geometry/building_blocks/Bypass.hpp>

using namespace RedMA;

int main(int argc, char **argv)
{
    #ifdef HAVE_MPI
    MPI_Init (nullptr, nullptr);
    shp<Epetra_Comm> comm (new Epetra_MpiComm(MPI_COMM_WORLD));
    #else
    shp<Epetra_Comm> comm(new Epetra_SerialComm ());
    #endif

    DataContainer data;
    data.setDatafile("data");
    data.setVerbose(comm->MyPID() == 0);

    const double A = data("parameters/stenosis_amplitude", 0.4);
    const double W = data("parameters/stenosis_width", 0.15);
    const unsigned int N = data("parameters/stenosis_number", 0);

    // Tube tb(comm, "fine", true, 1, 2);
    // BifurcationSymmetric tb(comm, "fine", true, 90);
    Bypass tb(comm, "bypass", true, true, true, N, true);

    tb.setDatafile(data);

    tb.readMesh();
    tb.setParameterValue("stenosis_amplitude", A);
    tb.setParameterValue("stenosis_width", W);
    tb.applyGlobalTransformation();

    tb.dumpMesh("output/", "../../../meshes/", "bypass");


    return 0;
}
