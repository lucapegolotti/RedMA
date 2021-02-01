#include "ProblemRB.hpp"

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
    typedef BlockAssembler BAssembler;
    printlog(MAGENTA, "[ProblemRB] starting setup ... \n", M_data.getVerbose());

    std::string geometriesDir = M_data("geometric_structure/geometries_dir",
                                       "../../../meshes/");

    M_tree.readMeshes(geometriesDir);
    M_tree.traverseAndDeformGeometries();

    M_assembler.reset(new BAssembler(M_data, M_tree));
    M_TMAlgorithm = TimeMarchingAlgorithmFactory(M_data, M_assembler);
    printlog(MAGENTA, "done\n", M_data.getVerbose());
}

void
ProblemRB::
solve()
{
    double t0 = M_data("time_discretization/t0", 0.0);
    double t0ramp = M_data("time_discretization/t0ramp", 0.0);
    double T = M_data("time_discretization/T", 1.0);
    double dt = M_data("time_discretization/dt", 0.01);
    int saveEvery = M_data("exporter/save_every", 1);

    double t;

    if (t0 - t0ramp > 1e-15)
        t = t0ramp;
    else
        t = t0;

    unsigned int count = 0;

    while (T - t > dt/2)
    {
        if (t < t0)
        {
            std::string msg = "[ProblemRB] solving ramp"
                              ", t = " + std::to_string(t) + " -> " + std::to_string(t+dt) + "\n";
            printlog(MAGENTA, msg, M_data.getVerbose());
        }
        else
        {
            std::string msg = "[ProblemRB] solving timestep " + std::to_string(count) +
                              ", t = " + std::to_string(t) + " -> " + std::to_string(t+dt) + "\n";
            printlog(MAGENTA, msg, M_data.getVerbose());
        }

        int status = -1;

        M_solution = spcast<BlockVector>(M_TMAlgorithm->advance(t, dt, status));

        if (status)
            throw new Exception("Error in solver. Status != 0.");

        t += dt;

        if (M_storeSolutions && t >= t0)
        {
            M_solutions.push_back(M_solution);
            M_timestepsSolutions.push_back(t);
        }

        if (M_storeNonLinearTerms && t >= t0)
        {
            BBV nonLinearTerm = spcast<BlockVector>(M_assembler->getNonLinearTerm());
            M_nonLinearTerms.push_back(nonLinearTerm);
        }

        if (t >= t0 && saveEvery > 0 && count % saveEvery == 0)
            M_assembler->exportSolution(t, M_solution);
        M_assembler->postProcess(t, dt, M_solution);
        M_TMAlgorithm->shiftSolutions(M_solution);
        if (t >= t0)
            count++;
    }
}

}
