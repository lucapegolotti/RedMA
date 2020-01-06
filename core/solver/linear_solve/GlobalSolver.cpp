#include <GlobalSolver.hpp>

namespace RedMA
{

GlobalSolver::
GlobalSolver(const GetPot& datafile, commPtr_Type comm, bool verbose,
             AbstractFunctor* exactSolution) :
  M_geometryParser(datafile("geometric_structure/xmlfile","tree.xml"),
                   comm, verbose),
  M_datafile(datafile),
  M_comm(comm),
  M_globalAssembler(M_datafile, M_comm, M_verbose),
  M_exportNorms(false),
  M_exportErrors(false),
  M_verbose(verbose)
{
    M_tree = M_geometryParser.getTree();

    std::string geometriesDir = datafile("geometric_structure/geometries_dir",
                                         "../../../meshes/");

    M_tree.readMeshes(geometriesDir);
    M_tree.traverseAndDeformGeometries();

    M_globalAssembler.setTreeStructure(M_tree);
    M_globalAssembler.setup();

    // we set it here because boundary conditions might be set in the constructor
    // of the time advancing scheme
    if (exactSolution != nullptr)
        setExactSolution(exactSolution);

    M_timeMarchingAlgorithm =
            TimeMarchingAlgorithmsFactory(datafile,
                                                         &M_globalAssembler,
                                                         M_comm, M_verbose);

    M_globalAssembler.setTimeIntegrationOrder(M_timeMarchingAlgorithm->getOrder());
}

void
GlobalSolver::
solve()
{
    std::ofstream outFile;

    if (M_exportNorms)
    {
        outFile.open(M_filename);
        // outFile << AbstractAssembler::normFileFirstLine() << std::flush;
    }
    if (M_exportErrors)
    {
        outFile.open(M_filename);
        // outFile << AbstractAssembler::errorsFileFirstLine() << std::flush;
    }

    bool steady = M_datafile("time_discretization/steady", false);
    if (steady)
    {
        TimeMarchingAlgorithmPtr hdlrAlgorithm = M_timeMarchingAlgorithm;

        double t0 = 0.0;
        double dt = 0.0;

        // hdlrAlgorithm->setInitialCondition(M_globalAssembler.getInitialCondition());

        solveTimestep(t0, dt);

        M_globalAssembler.setTimeAndPrevSolution(t0,
                                                 hdlrAlgorithm->getSolution());
        M_globalAssembler.postProcess();

        M_globalAssembler.exportSolutions(t0, hdlrAlgorithm->getSolution());

        if (M_exportNorms)
        {
            M_globalAssembler.appendNormsToFile(t0, hdlrAlgorithm->getSolution(),
                                                outFile);
        }
        if (M_exportErrors)
        {
            M_globalAssembler.appendErrorsToFile(t0, hdlrAlgorithm->getSolution(),
                                                 outFile);
        }

    }
    else
    {
        double t0 = M_datafile("time_discretization/t0", 0.0);
        double T = M_datafile("time_discretization/T", 1.0);
        double dt = M_datafile("time_discretization/dt", 0.01);
        unsigned int save_every = M_datafile("exporter/save_every", 1);

        double t = t0;
        TimeMarchingAlgorithmPtr hdlrAlgorithm = M_timeMarchingAlgorithm;
        hdlrAlgorithm->setInitialCondition(M_globalAssembler.getInitialCondition());
        M_globalAssembler.exportSolutions(t, hdlrAlgorithm->getSolution());
        if (M_exportErrors)
        {
            M_globalAssembler.appendErrorsToFile(t, hdlrAlgorithm->getSolution(),
                                                 outFile);
        }

        unsigned int count = 1;
        while (T - t > dt/2)
        {
            M_globalAssembler.setTimestep(dt);
            solveTimestep(t, dt);
            t += dt;
            // we do this so that the single assemblers are aware of the
            // solution (e.g. for post processing)
            M_globalAssembler.setTimeAndPrevSolution(t,
                                                     hdlrAlgorithm->getSolution());
            M_globalAssembler.postProcess();
            if (count % save_every == 0)
                M_globalAssembler.exportSolutions(t, hdlrAlgorithm->getSolution());
            if (M_exportNorms)
            {
                M_globalAssembler.appendNormsToFile(t, hdlrAlgorithm->getSolution(),
                                                    outFile);
            }
            if (M_exportErrors)
            {
                M_globalAssembler.appendErrorsToFile(t, hdlrAlgorithm->getSolution(),
                                                     outFile);
            }
            count++;
        }
    }
    if (M_exportNorms || M_exportErrors)
        outFile.close();
}

void
GlobalSolver::
setExportNorms(std::string filename)
{
    M_exportNorms = true;
    M_filename = filename;
}

void
GlobalSolver::
setExportErrors(std::string filename)
{
    M_exportErrors = true;
    M_filename = filename;
}

void
GlobalSolver::
setForcingFunction(FunctionType forcingFunction,
                   FunctionType forcingFunctionDt)
{
    M_globalAssembler.setForcingFunction(forcingFunction, forcingFunctionDt);
}

void
GlobalSolver::
solveTimestep(const double& time, double& dt)
{
    M_timeMarchingAlgorithm->solveTimestep(time, dt);
}

void
GlobalSolver::
setLawInflow(std::function<double(double)> maxLaw)
{
    M_globalAssembler.setLawInflow(maxLaw);
}

void
GlobalSolver::
setExactSolution(AbstractFunctor* exactSolution)
{
    M_globalAssembler.setExactSolution(exactSolution);
}

void
GlobalSolver::
setLawDtInflow(std::function<double(double)> maxLawDt)
{
    M_globalAssembler.setLawDtInflow(maxLawDt);
}

void
GlobalSolver::
printMeshSize(std::string filename)
{
    M_globalAssembler.printMeshSize(filename);
}

}  // namespace RedMA
