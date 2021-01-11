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

#include <redma/geometry/Tube.hpp>
#include <redma/geometry/BifurcationSymmetric.hpp>

using namespace RedMA;

int main(int argc, char **argv)
{
    #ifdef HAVE_MPI
    MPI_Init (nullptr, nullptr);
    shp<Epetra_Comm> comm (new Epetra_MpiComm(MPI_COMM_WORLD));
    #else
    shp<Epetra_Comm> comm(new Epetra_SerialComm ());
    #endif

    Tube tb(comm, "fine", true, 1, 2);
    // BifurcationSymmetric tb(comm, "fine", true, 90);

    // M_parametersHandler.registerParameter("bend", 0.0, 0, M_PI/2, randomizible, false);
    // M_parametersHandler.registerParameter("L_ratio", 1.0, 0.7, 1.3, randomizible);
    // M_parametersHandler.registerParameter("Rout_ratio", 1.0, 0.6, 1.4, randomizible);

    tb.readMesh();
    tb.setParameterValue("Rout_ratio", 0.6);
    tb.setParameterValue("bend", M_PI*0.3);
    tb.applyGlobalTransformation();
    tb.dumpMesh("output/", "../../../meshes/", "tube");


    return 0;
}
