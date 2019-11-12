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
#include <Tube.hpp>

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
    bifurcation.setParameterValue("alpha", 2.7665588261156668);
    bifurcation.setParameterValue("alpha_axis", -1.7026886923864399);
    bifurcation.setParameterValue("bx", 0.7042145287300223);
    bifurcation.setParameterValue("by", 3.143530101569969);
    bifurcation.setParameterValue("bz", -7.4118393744460951);
    double a = 0.7;
    double b = 0.3;
    double c = std::sqrt(1.0 - a*a - b*b);
    bifurcation.setParameterValue("rotation_axis_x", a);
    bifurcation.setParameterValue("rotation_axis_y", b);
    bifurcation.setParameterValue("rotation_axis_z", c);
    bifurcation.setParameterValue("scale", 1);

    bifurcation.setParameterValue("out1_alphax", -2);
    bifurcation.setParameterValue("out1_alphay", -2);
    bifurcation.setParameterValue("out1_alphaz", -2);
    bifurcation.setParameterValue("out2_alphax", -2);
    bifurcation.setParameterValue("out2_alphay", -2);
    bifurcation.setParameterValue("out2_alphaz", -2);

    bifurcation.applyGlobalTransformation();
    bifurcation.dumpMesh("output/", "../../../meshes/", "deformedBifurcation");
    return 0;
}
