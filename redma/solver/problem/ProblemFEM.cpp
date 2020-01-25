#include "redma/solver/problem/ProblemFEM.hpp"

namespace RedMA
{

ProblemFEM::
ProblemFEM(const GetPot& datafile, EPETRACOMM comm, bool verbose) :
  aProblem(datafile),
  M_geometryParser(datafile,
                   datafile("geometric_structure/xmlfile","tree.xml"),
                   comm, verbose)
{
    M_tree = M_geometryParser.getTree();

    std::string geometriesDir = datafile("geometric_structure/geometries_dir",
                                         "../../../meshes/");

    M_tree.readMeshes(geometriesDir);
    M_tree.traverseAndDeformGeometries();

    setup();
}

void
ProblemFEM::
setup()
{
    typedef BlockAssembler<BV, BM> BAssembler;
    M_TMAlgorithm = TimeMarchingAlgorithmFactory<BV, BM>(M_datafile);
    M_assembler.reset(new BAssembler(M_datafile, M_tree));
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
    int status = -1;
    M_solution = M_TMAlgorithm->advance(t, dt, M_assembler, status);
}

}
