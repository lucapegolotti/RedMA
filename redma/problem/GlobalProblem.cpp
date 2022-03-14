#include "GlobalProblem.hpp"

namespace RedMA
{

GlobalProblem::
GlobalProblem(const DataContainer& data, EPETRACOMM comm, bool doSetup) :
  aProblem(data),
  M_geometryParser(data,
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

    if (M_data("time_discretization/order", 0) > 0)
    {
        // initialize time marching algorithm (for instance BDF)
        M_TMAlgorithm = TimeMarchingAlgorithmFactory(M_data, M_assembler);
        M_TMAlgorithm->setComm(M_comm);
    }
    else
    {
        // initialize solver for the steady problem
        M_steadySolver.reset(new SteadySolver(M_data, M_assembler));
        M_steadySolver->setComm(M_comm);

        std::string IC_path = M_data("newton_method/ic_path", "");
        M_steadySolver->setInitialGuess(IC_path);
        shp<aVector> tmp = M_steadySolver->getInitialGuess();
        M_assembler->applyPiola(tmp, false);
        M_steadySolver->setInitialGuess(tmp);
    }

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

void
GlobalProblem::
solveSteady()
{
    std::string msg = "[GlobalProblem] solving steady problem...\n";
    printlog(MAGENTA, msg, M_data.getVerbose());

    int status = -1;

    M_solution = spcast<BlockVector>(M_steadySolver->solve(status));
    if (status)
        throw new Exception("Error in solver. Status != 0");

    M_assembler->exportSolution(0.0, M_solution);
    if (M_data("exporter/exporttotxt", true))
    {
        std::string fname = M_data("exporter/pathtotxt", "IC/");
        M_assembler->exportSolutionToTxt(0.0, M_solution, fname);
    }

}

void
GlobalProblem::
exportFromFiles(const std::string &inPath)
{
    double t0 = M_data("time_discretization/t0", 0.0);
    double T = M_data("time_discretization/T", 1.0);
    double dt = M_data("time_discretization/dt", 0.01);
    int saveEvery = M_data("exporter/save_every", 1);

    double t = t0;
    unsigned int count = 0;

    std::map<unsigned int, std::vector<shp<aVector>>> solutions = M_assembler->importSolution(inPath);

    while (T - t > dt/2)
    {
        std::string msg = "[GlobalProblem] exporting solution at timestep " + std::to_string(count) +
                          ", t = " + std::to_string(t+dt) + "\n";
        printlog(MAGENTA, msg, M_data.getVerbose());

        M_solution = spcast<BlockVector>(M_assembler->getZeroVector());

        unsigned int nBlocks = solutions.size();
        for (unsigned int innerCount=0; innerCount<nBlocks; ++innerCount)
        {
            spcast<BlockVector>(M_solution->block(innerCount))->setBlock(0,
                                                                         spcast<DistributedVector>(spcast<BlockVector>(solutions[innerCount][count])->block(0)));
            spcast<BlockVector>(M_solution->block(innerCount))->setBlock(1,
                                                                         spcast<DistributedVector>(spcast<BlockVector>(solutions[innerCount][count])->block(1)));
        }

        t += dt;

        if (t >= t0 && saveEvery > 0 && count % saveEvery == 0)
            M_assembler->exportSolution(t, M_solution);

        if (t >= t0)
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
