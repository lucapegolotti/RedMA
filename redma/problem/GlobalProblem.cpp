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
    int N_periods = M_data("time_discretization/N_periods", 1);
    double dt = M_data("time_discretization/dt", 0.01);
    int saveEvery = M_data("exporter/save_every", 1);
    int saveRamp = M_data("exporter/save_ramp", 1);

    double t;

    if (t0 - t0ramp > 1e-15)
        t = t0ramp;
    else
        t = t0;

    unsigned int count = 1;
    while (T * N_periods - t > dt/2)
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

        if (M_storeSolutions)
        {
            if (saveRamp || (t > t0))
            {
                M_solutions.push_back(M_solution);
                M_timestepsSolutions.push_back(t);
            }
            if ((std::abs(t-t0) < dt/2) || (std::abs(t-(t0-dt)) < dt/2))
                M_initialConditions.push_back(M_solution);
        }

        M_assembler->postProcess(t, M_solution);

        // saving displacement field, if necessary
        if (!(std::strcmp(M_data("assembler/type", "navierstokes").c_str(), "navierstokes_membrane")))
        {
            if (M_storeSolutions)
            {
                shp<BlockVector> tmpDisplacement (new BlockVector(0));
                tmpDisplacement->deepCopy(M_assembler->getDisplacement());
                if (saveRamp || (t > t0))
                    M_extraSolutions.push_back(tmpDisplacement);
                if ((std::abs(t-t0) < dt/2) || (std::abs(t-(t0-dt)) < dt/2))
                    M_extraInitialConditions.push_back(spcast<BlockVector>(M_assembler->getDisplacement()));
            }
        }

        if ((t > t0 && saveEvery > 0 && count % saveEvery == 0) || (std::abs(t-t0) < dt/2))
            M_assembler->exportSolution(t, M_solution);

        M_TMAlgorithm->shiftSolutions(M_solution);

        if ((t-t0) > dt/2)
            count++;
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

    std::map<unsigned int, std::vector<shp<BlockVector>>> solutions = M_assembler->importSolution(inPath);

    while (T - t > dt/2)
    {
        std::string msg = "[GlobalProblem] exporting solution at timestep " + std::to_string(count) +
                          ", t = " + std::to_string(t+dt) + "\n";
        printlog(MAGENTA, msg, M_data.getVerbose());

        unsigned int nBlocks = solutions.size();
        unsigned int nFields = solutions[0][0]->nRows();

        M_solution.reset(new BlockVector(nBlocks));

        for (unsigned int innerCount=0; innerCount<nBlocks; ++innerCount)
        {
            shp<BlockVector> tmpBlock(new BlockVector(nFields));
            for (unsigned int k=0; k<nFields; ++k)
                tmpBlock->setBlock(k, spcast<DistributedVector>(solutions[innerCount][count]->block(k)));

            M_solution->setBlock(innerCount, tmpBlock);
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

void
GlobalProblem::
validateRBConvectiveTerm()
{
    std::string filename = "ConvectiveTerm/RB_solution";
    M_assembler->getAssemblersMap()[0]->computeConvectiveTermFromFile(filename);
}
}
