#include "redma/solver/problem/ProblemFEM.hpp"

namespace RedMA
{

ProblemFEM::
ProblemFEM(const GetPot& datafile) :
  aProblem(datafile)
{
    setup();
}

void
ProblemFEM::
setup()
{
    M_timeMarchingAlgorithm =
                    TimeMarchingAlgorithmFactory<FEVECTOR,FEMATRIX>(M_datafile);
    M_assembler = AssemblerFactory<FEVECTOR,FEMATRIX>(M_datafile);
}

void
ProblemFEM::
solve()
{
    double t0 = M_datafile("time_discretization/t0", 0.0);
    double T = M_datafile("time_discretization/T", 1.0);
    double dt = M_datafile("time_discretization/dt", 0.01);
    unsigned int saveEvery = M_datafile("exporter/save_every", 1);

    double t = t;

    unsigned int count = 1;
    while (T - t > t/2)
    {
        solveTimestep(t, dt);
        t += dt;

        if (count % saveEvery == 0)
            M_assembler->exportSolution(t);

        M_assembler->postProcess();
        count++;
    }
}

void
ProblemFEM::
solveTimestep(const double& t, double& dt)
{
    M_solution = M_timeMarchingAlgorithm->advance(t, dt, M_assembler);
}

}
