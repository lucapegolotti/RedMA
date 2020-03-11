#include "redma/solver/problem/ProblemRB.hpp"

namespace RedMA
{

ProblemRB::
ProblemRB(const DataContainer& data, EPETRACOMM comm, bool doSetup) :
  aProblem(data),
  M_geometryParser(data.getDatafile(),
                   data("geometric_structure/xmlfile","tree.xml"),
                   comm, data.getVerbose()),
  M_storeSolutions(false)
{
    M_tree = M_geometryParser.getTree();

    if (doSetup)
        setup();
}

void
ProblemRB::
setup()
{
    typedef BlockAssembler<BV, BM> BAssembler;
    printlog(MAGENTA, "[ProblemRB] starting setup ... \n", M_data.getVerbose());

    std::string geometriesDir = M_data("geometric_structure/geometries_dir",
                                       "../../../meshes/");

    M_tree.readMeshes(geometriesDir);
    M_tree.traverseAndDeformGeometries();

    M_assembler.reset(new BAssembler(M_data, M_tree));
    M_TMAlgorithm = TimeMarchingAlgorithmFactory<BV, BM>(M_data, M_assembler);
    printlog(MAGENTA, "done\n", M_data.getVerbose());
}

void
ProblemRB::
solve()
{
    double t0 = M_data("time_discretization/t0", 0.0);
    double T = M_data("time_discretization/T", 1.0);
    double dt = M_data("time_discretization/dt", 0.01);
    unsigned int saveEvery = M_data("exporter/save_every", 1);

    double t = t0;

    unsigned int count = 1;
    while (T - t > dt/2)
    {
        std::string msg = "[ProblemFEM] solving timestep " + std::to_string(count) +
                          ", t = " + std::to_string(t) + " -> " + std::to_string(t+dt) + "\n";
        printlog(MAGENTA, msg, M_data.getVerbose());
        int status = solveTimestep(t, dt);
        if (status)
            throw new Exception("Error in solver. Status != 0.");

        t += dt;

        if (count % saveEvery == 0)
            M_assembler->exportSolution(t, M_solution);
        M_assembler->postProcess(t, M_solution);
        M_TMAlgorithm->shiftSolutions(M_solution);
        count++;
    }
}

int
ProblemRB::
solveTimestep(const double& t, double& dt)
{
    int status = -1;
    M_solution = M_TMAlgorithm->advance(t, dt, status);
    if (M_storeSolutions)
        M_solutions.push_back(M_solution);
    return status;
}

}
