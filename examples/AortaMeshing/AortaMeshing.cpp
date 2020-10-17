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

#include <redma/geometry/SegmentationsMerger.hpp>
#include <redma/geometry/GeometryPrinter.hpp>

#include <lifev/core/filter/GetPot.hpp>

using namespace RedMA;

int main(int argc, char **argv)
{
    #ifdef HAVE_MPI
    MPI_Init (nullptr, nullptr);
    std::shared_ptr<Epetra_Comm> comm (new Epetra_MpiComm(MPI_COMM_WORLD));
    #else
    std::shared_ptr<Epetra_Comm> comm(new Epetra_SerialComm ());
    #endif

    GetPot datafile("datafiles/data");

    SegmentationsMerger merger(datafile,comm,true);

    std::shared_ptr<SegmentationParser> sp1(new SegmentationParser(datafile, comm,
                "datafiles/aorta.pth", "datafiles/aorta_segs_final.ctgr", "linear", true));

    double constCenters = datafile("segmentation_parser/const_centers", 2.0);
    double constNormal = datafile("segmentation_parser/const_normals", 1.0);

    auto contours = sp1->getContours();

    std::cout << "contours = " << contours.size() << std::endl << std::flush;

    GeometricFace final;

    contours[std::atoi(argv[1])].print();


    TreeStructure tree1 = sp1->createTreeForward(1, constCenters, constNormal,
                                                nullptr, &contours[80]);

    GeometryPrinter printer;
    printer.saveToFile(tree1, "tree_1.xml", comm);

    TreeStructure tree2 = sp1->createTreeForward(3, constCenters, constNormal,
                                                &contours[111], nullptr);
    printer.saveToFile(tree2, "tree_2.xml", comm);
    tree2.readMeshes("../../../meshes/");
    tree2.traverseAndDeformGeometries();
    tree2.dump("output/","../../../meshes/");

    return 0;
}
