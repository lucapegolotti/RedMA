#include "GlobalProblem.hpp"

namespace RedMA
{

GlobalProblem::
GlobalProblem(const DataContainer& data, EPETRACOMM comm, bool doSetup) :
  aProblem(data),
  M_geometryParser(data.getDatafile(),
                   data("geometric_structure/xmlfile","tree.xml"),
                   comm, data.getVerbose()),
  M_storeSolutions(false),
  M_comm(comm)
{
    M_tree = M_geometryParser.getTree();

    if (doSetup)
        setup();
}

void
GlobalProblem::
setup()
{
    typedef BlockAssembler BAssembler;
    printlog(MAGENTA, "[Problem] starting setup ... \n", M_data.getVerbose());

    std::string geometriesDir = M_data("geometric_structure/geometries_dir",
                                       "../../../meshes/");

    // read and deform meshes contained in the tree.xml file
    M_tree.readMeshes(geometriesDir);
    M_tree.traverseAndDeformGeometries();

    // uncomment this to dump tree after it has been read
    // M_tree.dump("tree/","../../../meshes/");

    // create default assemblers. These are needed when we require the map from
    // the reference configuration to the deformed one (e.g., for Piola)
    M_defaultAssemblers.reset(new DefaultAssemblersLibrary
                             (M_data, M_tree.getMeshListNames(), M_comm));

    // create block assembler
    M_assembler.reset(new BAssembler(M_data, M_tree, M_defaultAssemblers));

    // initialize time marching algorithm (for instance BDF)
    M_TMAlgorithm = TimeMarchingAlgorithmFactory(M_data, M_assembler);
    M_TMAlgorithm->setComm(M_comm);
    printlog(MAGENTA, "done\n", M_data.getVerbose());
}

void
GlobalProblem::
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

    unsigned int count = 1;
    while (T - t > dt/2)
    {
        if (t < t0)
        {
            std::string msg = "[GlobalProblem] solving ramp"
                              ", t = " + std::to_string(t) + " -> " + std::to_string(t+dt) + "\n";
            printlog(MAGENTA, msg, M_data.getVerbose());
        }
        else
        {
            std::string msg = "[GlobalProblem] solving timestep " + std::to_string(count) +
                              ", t = " + std::to_string(t) + " -> " + std::to_string(t+dt) + "\n";
            printlog(MAGENTA, msg, M_data.getVerbose());
        }

        int status = -1;

        M_solution = spcast<BlockVector>(M_TMAlgorithm->advance(t, dt, status));
        if (status)
            throw new Exception("Error in solver. Status != 0");

        t += dt;

        if (M_storeSolutions && t >= t0)
        {
            M_solutions.push_back(M_solution);
            M_timestepsSolutions.push_back(t);
        }

        M_assembler->postProcess(t, M_solution);

        if ((t > t0 && saveEvery > 0 && count % saveEvery == 0) || (std::abs(t-t0) < dt/2))
            M_assembler->exportSolution(t, M_solution);

        M_TMAlgorithm->shiftSolutions(M_solution);

        if ((t-t0) > dt/2)
            count++;
    }
}

bool
GlobalProblem::
isFEProblem()
{
    return M_assembler->arePrimalAssemblersFE();
}

bool
GlobalProblem::
isRBProblem()
{
    return M_assembler->arePrimalAssemblersRB();
}

}
