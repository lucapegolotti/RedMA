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
    shp<Epetra_Comm> comm (new Epetra_MpiComm(MPI_COMM_WORLD));
    #else
    shp<Epetra_Comm> comm(new Epetra_SerialComm ());
    #endif

    GetPot datafile("datafiles/data");

    SegmentationsMerger merger(datafile,comm,true);

    shp<SegmentationParser> sp1(new SegmentationParser(datafile, comm,
                "datafiles/aorta.pth", "datafiles/aorta_segs_final.ctgr", "linear", true));

    double constCenters = datafile("segmentation_parser/const_centers", 2.0);
    double constNormal = datafile("segmentation_parser/const_normals", 1.0);

    auto contours = sp1->getContours();

    std::cout << "contours = " << contours.size() << std::endl << std::flush;

    std::vector<GeometricFace> cuts;
    // contours[std::atoi(argv[1])].print();


    TreeStructure tree1 = sp1->createTreeForward(1, constCenters, constNormal,
                                                nullptr, &contours[80]);

    cuts.push_back(contours[80]);

    GeometryPrinter printer;
    printer.saveToFile(tree1, "tree_aorta1.xml", comm);

    TreeStructure tree2 = sp1->createTreeForward(3, constCenters, constNormal,
                                                &contours[111], nullptr);

    cuts.push_back(contours[111]);

    printer.saveToFile(tree2, "tree_aorta2.xml", comm);

    shp<SegmentationParser> sp2(new SegmentationParser(datafile, comm,
                "datafiles/subclavian.pth", "datafiles/subclavian_segs_final.ctgr", "linear", true));

    contours = sp2->getContours();
    TreeStructure tree_subclavian = sp2->createTreeForward(2, constCenters, constNormal,
                                                           &contours[30], nullptr);

    cuts.push_back(contours[30]);

    printer.saveToFile(tree_subclavian, "tree_subclavian.xml", comm);
    // tree_subclavian.readMeshes("../../../meshes/");
    // tree_subclavian.traverseAndDeformGeometries();
    // tree_subclavian.dump("output/","../../../meshes/");

    shp<SegmentationParser> sp3(new SegmentationParser(datafile, comm,
                "datafiles/Lt-carotid.pth", "datafiles/Lt-carotid_segs_final.ctgr", "linear", true));

    contours = sp3->getContours();
    TreeStructure tree_ltcarotid = sp3->createTreeForward(2, constCenters, constNormal,
                                                           &contours[26], nullptr);

    cuts.push_back(contours[26]);

    shp<SegmentationParser> sp4(new SegmentationParser(datafile, comm,
               "datafiles/btrunk.pth", "datafiles/btrunk_segs_final.ctgr", "linear", true));

    contours = sp4->getContours();
    TreeStructure tree_btrunk = sp4->createTreeForward(2, constCenters, constNormal,
                                                          &contours[26], &contours[64]);

    cuts.push_back(contours[26]);
    cuts.push_back(contours[64]);

    printer.saveToFile(tree_btrunk, "tree_btrunk.xml", comm);
    // tree_btrunk1.readMeshes("../../../meshes/");
    // tree_btrunk1.traverseAndDeformGeometries();
    // tree_btrunk1.dump("output/","../../../meshes/");

    shp<SegmentationParser> sp5(new SegmentationParser(datafile, comm,
               "datafiles/rt-carotid.pth", "datafiles/rt-carotid_segs_final.ctgr", "linear", true));

    contours = sp5->getContours();
    TreeStructure tree_rtcarotid = sp5->createTreeForward(2, constCenters, constNormal,
                                                          &contours[80], nullptr);

    cuts.push_back(contours[80]);

    printer.saveToFile(tree_rtcarotid, "tree_rtcarotid.xml", comm);
    tree_rtcarotid.readMeshes("../../../meshes/");
    tree_rtcarotid.traverseAndDeformGeometries();
    tree_rtcarotid.dump("output/","../../../meshes/");

    for (auto cut : cuts)
        cut.print();

    return 0;
}
