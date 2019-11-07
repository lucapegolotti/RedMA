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

#include <tinyxml2.h>
#include <vector>

#include <GeometryPrinter.hpp>
#include <GeometryParser.hpp>
#include <BifurcationSymmetric.hpp>

using namespace RedMA;

int main(int argc, char **argv)
{
    #ifdef HAVE_MPI
    MPI_Init (nullptr, nullptr);
    std::shared_ptr<Epetra_Comm> comm (new Epetra_MpiComm(MPI_COMM_WORLD));
    #else
    std::shared_ptr<Epetra_Comm> comm(new Epetra_SerialComm ());
    #endif

    BifurcationSymmetric bifurcation(comm, "fine", true);
    bifurcation.readMesh("../../../meshes/");
    bifurcation.setParameterValue("out1_alphax", 0.0);
    bifurcation.setParameterValue("out1_alphaz", 0.5);
    bifurcation.setParameterValue("out2_alphaz", 0.5);
    //bifurcation.setParameterValue("out2_alphaz", 0.5);
    bifurcation.setParameterValue("out2_alphax", 0.5);
    // bifurcation.setParameterValue("out2_alphax", 0.2);
    // bifurcation.setParameterValue("out2_alpha_plane", -0.1);
    bifurcation.applyGlobalTransformation();
    bifurcation.dumpMesh("output/", "../../../meshes/", "deformedBifurcation");
    return 0;
}
