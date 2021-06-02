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
            throw new Exception("Error in solver. Status != 0.");

        t += dt;

        if (M_storeSolutions && t >= t0)
        {
            M_solutions.push_back(M_solution);
            M_timestepsSolutions.push_back(t);
        }
        if (t >= t0 && saveEvery > 0 && count % saveEvery == 0)
            M_assembler->exportSolution(t, M_solution);

        M_assembler->postProcess(t, M_solution);

        /*std::cout << M_assembler->getDualAssemblers()[0]->getChildBT()->multiplyByVector(M_solution->block(1))->norm2() << std::endl <<std::flush;
        std::cout << M_solution->block(1)->norm2() << std::endl << std::flush;

        shp<BlockMatrix> Mat1 = spcast<BlockMatrix>(M_assembler->getAssemblersMap()[0]->getMass(0.0, nullptr));
        shp<BlockMatrix> Mat2 = spcast<BlockMatrix>(M_assembler->getAssemblersMap()[0]->assembleMatrix(1));
        Mat2->multiplyByScalar(dt * 2.0 / 3.0);

        shp<BlockMatrix> Mat3 = spcast<BlockMatrix>(M_assembler->getAssemblersMap()[0]->assembleMatrix(2));
        Mat3->multiplyByScalar(dt * 2.0 / 3.0);

        Mat1->add(Mat2);
        Mat1->add(Mat3);

        shp<BlockMatrix> Mat4 = spcast<BlockMatrix>(M_assembler->getDualAssemblers()[0]->getChildBT());
        Mat4->multiplyByScalar(dt * 2.0 / 3.0);

        shp<BlockMatrix> Mat(new BlockMatrix(2,2));
        Mat->setBlock(0,0, Mat1);
        Mat->setBlock(0,1, Mat4);

        std::cout << Mat->multiplyByVector(M_solution)->norm2() <<std::endl << std::flush;*/

        exit(1);

        M_TMAlgorithm->shiftSolutions(M_solution);

        if (t >= t0)
            count++;
    }
}

void
GlobalProblem::
exportFromFiles(const std::string &path)
{
    double t0 = M_data("time_discretization/t0", 0.0);
    double T = M_data("time_discretization/T", 1.0);
    double dt = M_data("time_discretization/dt", 0.01);
    int saveEvery = M_data("exporter/save_every", 1);

    double t = t0;
    unsigned int count = 0;

    std::map<unsigned int, std::vector<std::pair<shp<VECTOREPETRA>, shp<VECTOREPETRA>>>> solutions = M_assembler->importSolution(path);

    while (T - t > dt/2)
    {
        std::string msg = "[GlobalProblem] exporting solution at timestep " + std::to_string(count) +
                          ", t = " + std::to_string(t+dt) + "\n";
        printlog(MAGENTA, msg, M_data.getVerbose());

        M_solution = spcast<BlockVector>(M_assembler->getZeroVector());

        unsigned int nBlocks = solutions.size();
        for (unsigned int innerCount=0; innerCount<nBlocks; ++innerCount)
        {
            (M_solution->block(innerCount)->block(0)->data()).reset(new VECTOREPETRA(*solutions[innerCount][count].first));
            (M_solution->block(innerCount)->block(1)->data()).reset(new VECTOREPETRA(*solutions[innerCount][count].second));
        }

        t += dt;

        if (M_storeSolutions && t >= t0)
        {
            M_solutions.push_back(M_solution);
            M_timestepsSolutions.push_back(t);
        }
        if (t >= t0 && saveEvery > 0 && count % saveEvery == 0)
            M_assembler->exportSolution(t, M_solution);

        M_assembler->postProcess(t, M_solution);
        M_TMAlgorithm->shiftSolutions(M_solution);

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
