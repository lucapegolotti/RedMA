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

#include <SegmentationsMerger.hpp>
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

    unsigned int constVector = 1.0;
    unsigned int constNormal = 1.0;
    // unsigned int begin = 16;
    // if (argc == 4)
    // {
    //     constVector = std::atoi(argv[1]);
    //     constNormal = std::atoi(argv[2]);
    //     begin = std::atoi(argv[3]);
    // }
    //
    // SegmentationParser sp1(comm, "datafiles/aorta.pth", "datafiles/aorta.ctgr",
    //                       "linear", true);
    // TreeStructure tree1 = sp1.createTree(3, constVector, constNormal, begin);
    // GeometryPrinter printer1;
    // printer1.saveToFile(tree1, "tree_aorta.xml", comm);
    // tree1.readMeshes("../../../meshes/");
    // tree1.traverseAndDeformGeometries();
    // tree1.dump("output_aorta/","../../../meshes/");
    //
    // SegmentationParser sp2(comm, "datafiles/right_iliac.pth", "datafiles/right_iliac.ctgr",
    //                       "linear", true);
    // TreeStructure tree2 = sp2.createTree(3,constVector, constNormal);
    // GeometryPrinter printer2;
    // printer2.saveToFile(tree2, "tree_right_iliac.xml", comm);
    // tree2.readMeshes("../../../meshes/");
    // tree2.traverseAndDeformGeometries();
    // tree2.dump("output_right_iliac/","../../../meshes/");

    unsigned int begin = 16;
    if (argc == 4)
    {
        constVector = std::atoi(argv[1]);
        constNormal = std::atoi(argv[2]);
        begin = std::atoi(argv[3]);
    }

    SegmentationsMerger merger(comm);

    std::shared_ptr<SegmentationParser> sp1(new SegmentationParser(comm,
                "datafiles/aorta.pth", "datafiles/aorta.ctgr", "linear", true));

    std::shared_ptr<SegmentationParser> sp2(new SegmentationParser(comm,
                "datafiles/right_iliac.pth", "datafiles/right_iliac.ctgr", "linear", true));

    std::vector<std::shared_ptr<SegmentationParser> > parsers;
    parsers.push_back(sp1);
    parsers.push_back(sp2);

    TreeStructure tree = merger.merge(parsers);

    GeometryPrinter printer;
    printer.saveToFile(tree, "tree.xml", comm);
    tree.readMeshes("../../../meshes/");
    tree.traverseAndDeformGeometries();
    tree.dump("output/","../../../meshes/");

    // SegmentationParser sp1(comm, "datafiles/aorta.pth", "datafiles/aorta.ctgr",
    //                       "linear", true);
    // TreeStructure tree1 = sp1.createTree(3, constVector, constNormal, begin);
    // GeometryPrinter printer1;
    // printer1.saveToFile(tree1, "tree_aorta.xml", comm);
    // tree1.readMeshes("../../../meshes/");
    // tree1.traverseAndDeformGeometries();
    // tree1.dump("output_aorta/","../../../meshes/");
    //
    // SegmentationParser sp2(comm, "datafiles/right_iliac.pth", "datafiles/right_iliac.ctgr",
    //                       "linear", true);
    // TreeStructure tree2 = sp2.createTree(3,constVector, constNormal);
    // GeometryPrinter printer2;
    // printer2.saveToFile(tree2, "tree_right_iliac.xml", comm);
    // tree2.readMeshes("../../../meshes/");
    // tree2.traverseAndDeformGeometries();
    // tree2.dump("output_right_iliac/","../../../meshes/");


    return 0;
}
