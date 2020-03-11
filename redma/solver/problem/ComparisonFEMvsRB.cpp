#include "ComparisonFEMvsRB.hpp"

namespace RedMA
{

ComparisonFEMvsRB::
ComparisonFEMvsRB(const DataContainer& data, EPETRACOMM comm) :
  M_comm(comm),
  M_data(data),
  M_problemFEM(data, comm),
  M_problemRB(data,comm)
{

}

void
ComparisonFEMvsRB::
runFEM()
{
    LifeV::LifeChrono chrono;
    chrono.start();

    std::string msg = "Starting chrono\n";
    printlog(MAGENTA, msg, true);

    // deactivate exporter
    M_data.setValue("exporter/save_every", -1);

    M_problemFEM.doStoreSolutions();
    M_problemFEM.solve();

    M_timeFEM = chrono.diff();
    msg = "Total time =  ";
    msg += std::to_string(M_timeFEM);
    msg += " seconds\n";
    printlog(MAGENTA, msg, true);
}

void
ComparisonFEMvsRB::
runRB()
{
    LifeV::LifeChrono chrono;
    chrono.start();

    std::string msg = "Starting chrono\n";
    printlog(MAGENTA, msg, true);

    // deactivate exporter
    M_data.setValue("exporter/save_every", -1);

    M_problemRB.doStoreSolutions();
    M_problemRB.solve();

    M_timeRB = chrono.diff();
    msg = "Total time =  ";
    msg += std::to_string(M_timeRB);
    msg += " seconds\n";
    printlog(MAGENTA, msg, true);
}

void
ComparisonFEMvsRB::
postProcess()
{
    double t0 = M_data("time_discretization/t0", 0.0);
    double T = M_data("time_discretization/T", 1.0);
    double dt = M_data("time_discretization/dt", 0.01);

    double t = t0;

    unsigned int nsolutions = M_problemFEM.getSolutions().size();
    for (unsigned int i = 0; i < nsolutions; i++)
    {
        t = t + dt;
        auto diff = M_problemFEM.getSolutions()[i] -
        M_problemRB.getBlockAssembler()->convertFunctionRBtoFEM(M_problemRB.getSolutions()[i], M_comm);
        M_problemFEM.getBlockAssembler()->exportSolution(t, diff);
    }
}

}
