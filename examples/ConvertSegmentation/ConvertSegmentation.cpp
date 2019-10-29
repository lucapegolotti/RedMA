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

#include <SegmentationParser.hpp>
#include <GeometryPrinter.hpp>


using namespace RedMA;

int main(int argc, char **argv)
{
    #ifdef HAVE_MPI
    MPI_Init (nullptr, nullptr);
    std::shared_ptr<Epetra_Comm> comm (new Epetra_MpiComm(MPI_COMM_WORLD));
    #else
    std::shared_ptr<Epetra_Comm> comm(new Epetra_SerialComm ());
    #endif

    SegmentationParser sp(comm, "datafiles/aorta.pth", "datafiles/aorta.ctgr",
                          "linear", true);

    unsigned int begin = 16;
    if (argc == 2)
        begin = std::atoi(argv[1]);

    TreeStructure tree = sp.createTree(begin);
    GeometryPrinter printer;
    printer.saveToFile(tree, "tree.xml", comm);
    tree.readMeshes("../../../meshes/");
    tree.traverseAndDeformGeometries();
    tree.dump("output/","../../../meshes/");

    unsigned int count = 0;
    std::shared_ptr<TreeNode> p = tree.getRoot();
    while (p != nullptr)
    {
        std::cout << count << std::endl;
        p->M_block->getOutlet(0).print();
        p = p->M_children[0];
        count++;
    }

    return 0;
}
