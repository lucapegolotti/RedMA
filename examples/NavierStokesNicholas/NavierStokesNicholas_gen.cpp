#include <redma/RedMA.hpp>
#include <redma/problem/GlobalProblem.hpp>
#include <redma/problem/DataContainer.hpp>

#include <redma/reduced_basis/MatricesGenerator.hpp>

using namespace RedMA;

double inflow(double t, double a, double c)
{
    double T = 0.3;
    return 1-cos(2*M_PI*t/T) + c*sin(2*M_PI*a*t/T);
}

int main(int argc, char **argv)
{
    #ifdef HAVE_MPI
    MPI_Init (nullptr, nullptr);
    EPETRACOMM comm (new Epetra_MpiComm(MPI_COMM_WORLD));
    #else
    EPETRACOMM comm(new Epetra_SerialComm());
    #endif

    DataContainer data;
    data.setDatafile("datafiles/data_fem");
    data.setVerbose(comm->MyPID() == 0);

    MatricesGenerator generator(data, comm);
    generator.generate();

    return 0;
}
