//
// Created by micol on 26/04/2021.
//

#include <redma/RedMA.hpp>
#include <redma/problem/DataContainer.hpp>
#include <redma/reduced_basis/SnapshotsSampler.hpp>

using namespace RedMA;

double inflow2(double t)
{
    return 1;
    //return (t<0.02)*(1.0 - std::cos((t) * M_PI / (0.02 )) )/ 2.0 +(t>=0.02);
}

double inflow3(double t)
{
    return 1;
    //return (t<0.02)*(1.0 - std::cos((t) * M_PI / (0.02 )) )/ 2.0 +(t>=0.02);
}

int main(int argc, char **argv)
{
#ifdef HAVE_MPI
    MPI_Init(&argc, &argv);
    EPETRACOMM comm (new Epetra_MpiComm(MPI_COMM_WORLD));
#else
    EPETRACOMM comm(new Epetra_SerialComm());
#endif

    DataContainer data;
    data.setDatafile("datafiles/data");
    data.setVerbose(comm->MyPID() == 0);
    data.finalize();
    if (argc==3){
        double index=atof(argv[1]);
        double index_max=atof(argv[2]);
        SnapshotsSampler sampler(data, comm,0,inflow2,inflow3,index,index_max);
        sampler.takeSnapshots();
    }


    return 0;
}
