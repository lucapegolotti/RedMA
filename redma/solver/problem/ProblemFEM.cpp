#include "redma/solver/problem/ProblemFEM.hpp"

namespace RedMA
{

ProblemFEM::
ProblemFEM(const DataContainer& data, EPETRACOMM comm) :
  aProblem(data),
  M_geometryParser(data.getDatafile(),
                   data("geometric_structure/xmlfile","tree.xml"),
                   comm, data.getVerbose())
{
    M_tree = M_geometryParser.getTree();

    std::string geometriesDir = data("geometric_structure/geometries_dir",
                                     "../../../meshes/");

    M_tree.readMeshes(geometriesDir);
    M_tree.traverseAndDeformGeometries();

    setup();
}

ProblemFEM::
ProblemFEM(const DataContainer& data, const CommunicatorsDistributor& commD) :
  aProblem(data),
  M_geometryParser(data.getDatafile(),
                   data("geometric_structure/xmlfile","tree.xml"),
                   commD.getComms(), commD.getProcessMap(), data.getVerbose())
{
    M_tree = M_geometryParser.getTree();

    std::string geometriesDir = data("geometric_structure/geometries_dir",
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

    printlog(MAGENTA, "[ProblemFEM] starting setup ... \n", M_data.getVerbose());
    M_assembler.reset(new BAssembler(M_data, M_tree));
    M_TMAlgorithm = TimeMarchingAlgorithmFactory<BV, BM>(M_data, M_assembler);
    printlog(MAGENTA, "done\n", M_data.getVerbose());
}

void
ProblemFEM::
solve()
{
    double t0 = M_data("time_discretization/t0", 0.0);
    double T = M_data("time_discretization/T", 1.0);
    double dt = M_data("time_discretization/dt", 0.01);
    unsigned int saveEvery = M_data("exporter/save_every", 1);

    double t = t0;

    unsigned int count = 1;
    while (T - t > t/2)
    {
        std::string msg = "[ProblemFEM] solving timestep " + std::to_string(count-1) +
                          ", t = " + std::to_string(t) + "\n";
        printlog(MAGENTA, msg, M_data.getVerbose());
        solveTimestep(t, dt);
        t += dt;

        if (count % saveEvery == 0)
            M_assembler->exportSolution(t, M_solution);

        M_assembler->postProcess(M_solution);
        count++;
    }
}

void
ProblemFEM::
solveTimestep(const double& t, double& dt)
{
    int status = -1;
    M_solution = M_TMAlgorithm->advance(t, dt, status);
}

}
