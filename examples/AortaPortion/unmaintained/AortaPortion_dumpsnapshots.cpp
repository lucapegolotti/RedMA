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
#include <redma/reduced_basis/SnapshotsSampler.hpp>

using namespace RedMA;
// using namespace boost::filesystem;

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

    std::string snapshotsdir = data("rb/offline/snapshots/directory", "snapshots");

    if (!fs::exists(snapshotsdir))
        throw new Exception("Snapshots directory has not been generated yet!");

    std::string paramdir = snapshotsdir + "/param";
    unsigned int i = 0;
    // we loop over the folders with the parameters
    while (fs::exists(paramdir + std::to_string(i)))
    {
        if (fs::exists(paramdir + std::to_string(i) + "/tree.xml"))
        {
            GeometryParser gParser(data.getDatafile(),
                                   paramdir + std::to_string(i) + "/tree.xml", comm, true);
            comm->Barrier();

            TreeStructure& tree = gParser.getTree();
            tree.readMeshes("../../../meshes/");
            tree.traverseAndDeformGeometries();
            tree.dump(paramdir + std::to_string(i) + "/","../../../meshes/");
        }
        i++;
    }

    return 0;
}
