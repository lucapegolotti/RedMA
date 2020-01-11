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

#include <redma/geometry/GeometryPrinter.hpp>
#include <redma/geometry/GeometryParser.hpp>

using namespace RedMA;

int main(int argc, char **argv)
{
    #ifdef HAVE_MPI
    MPI_Init (nullptr, nullptr);
    std::shared_ptr<Epetra_Comm> comm (new Epetra_MpiComm(MPI_COMM_WORLD));
    #else
    std::shared_ptr<Epetra_Comm> comm(new Epetra_SerialComm ());
    #endif

    unsigned int nSub = 10;

    if (argc > 1)
        nSub = std::atoi(argv[1]);

    TreeStructure tree(true);
    tree.createRandom(nSub, comm);
    tree.readMeshes("../../../meshes/");
    tree.traverseAndDeformGeometries();
    tree.dump("output_original/","../../../meshes/");

    GeometryPrinter printer;
    printer.saveToFile(tree, "tree.xml", comm);

    GeometryParser gParser("tree.xml", comm, true);
    comm->Barrier();

    TreeStructure& tree2 = gParser.getTree();
    tree2.readMeshes("../../../meshes/");
    tree2.traverseAndDeformGeometries();
    tree2.dump("output_read/","../../../meshes/");

    return 0;
}
