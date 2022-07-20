//
// Created by Federico Betti on 10/04/2022.
//

#include <redma/RedMA.hpp>
#include <redma/problem/DataContainer.hpp>
#include <redma/reduced_basis/SnapshotsSteadySampler.hpp>

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
    // data.setValueString("geometric_structure/xmlfile", "tree_snapshots.xml");
    data.setVerbose(comm->MyPID() == 0);

    std::vector<unsigned int> numSamples = {1, 2, 1};

    unsigned int Nstart = 0;

    SnapshotsSteadySampler sampler(data, comm, numSamples);
    sampler.takeSnapshots(Nstart);

    return 0;
}