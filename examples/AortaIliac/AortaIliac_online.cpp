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

#include <redma/RedMA.hpp>
#include <redma/problem/DataContainer.hpp>
// #include <redma/problem/ProblemRB.hpp>

// #include <redma/problem/ComparisonFEMvsRB.hpp>

using namespace RedMA;

int main(int argc, char **argv)
{
    #ifdef HAVE_MPI
    MPI_Init (nullptr, nullptr);
    EPETRACOMM comm (new Epetra_MpiComm(MPI_COMM_WORLD));
    #else
    EPETRACOMM comm(new Epetra_SerialComm());
    #endif

    DataContainer data;
    data.setDatafile("datafiles/data");
    data.setVerbose(comm->MyPID() == 0);
    data.finalize();

    // GeometryParser gParser(data.getDatafile() , "geometries/tree.xml", comm, true);
    // comm->Barrier();
    //
    // TreeStructure& tree2 = gParser.getTree();
    // tree2.readMeshes("../../../meshes/");
    // tree2.traverseAndDeformGeometries();
    // tree2.dump("output_read/","../../../meshes/");

    // ComparisonFEMvsRB comparison(data, comm);
    //
    // std::string msg = "Starting chrono\n";
    // printlog(MAGENTA, msg, true);
    // Chrono chrono;
    // chrono.start();
    //
    // shp<ProblemRB> rbProblem(new ProblemRB(data, comm));
    // double setupTimeRB = chrono.diff();
    //
    // comparison.setProblemRB(rbProblem);
    //
    // comparison.runRB();
    // double runTimeRB = comparison.getTimeRB();
    // comparison.exportRB(4);
    //
    //
    // msg = "Total time =  ";
    // msg += std::to_string(chrono.diff());
    // msg += " seconds\n";
    // printlog(MAGENTA, msg, true);
    //
    // msg = "Setup time =  ";
    // msg += std::to_string(setupTimeRB);
    // msg += " seconds\n";
    // printlog(MAGENTA, msg, true);
    //
    // msg = "Run time =  ";
    // msg += std::to_string(runTimeRB);
    // msg += " seconds\n";
    // printlog(MAGENTA, msg, true);

    return 0;
}
