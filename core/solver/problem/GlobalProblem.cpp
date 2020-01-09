#include <GlobalProblem.hpp>

namespace RedMA
{

GlobalProblem::
GlobalProblem(const GetPot& datafile, commPtr_Type comm, bool verbose,
              AbstractFunctor* exactSolution) :
  M_datafile(datafile),
  M_comm(comm),
  M_exportNorms(false),
  M_exportErrors(false),
  M_verbose(verbose)
{
    // we set it here because boundary conditions might be set in the constructor
    // of the time advancing scheme
    if (exactSolution != nullptr)
        setExactSolution(exactSolution);

    setup();
}

void
GlobalProblem::
setup()
{
    M_assembler.reset(new GlobalAssembler(M_datafile, M_comm, M_verbose));
    std::string treeName = M_datafile("geometric_structure/xmlfile","tree.xml");
    M_geometryParser.reset(new GeometryParser(treeName, M_comm, M_verbose));

    M_tree = M_geometryParser->getTree();

    std::string geometriesDir = M_datafile("geometric_structure/geometries_dir",
                                           "../../../meshes/");

    M_tree.readMeshes(geometriesDir);
    M_tree.traverseAndDeformGeometries();

    M_assembler->setTreeStructure(M_tree);
    M_assembler->setup();

    M_timeMarchingAlgorithm =
            TimeMarchingAlgorithmsFactory(M_datafile, M_assembler.get(),
                                          M_comm, M_verbose);

    M_assembler->setTimeIntegrationOrder(M_timeMarchingAlgorithm->getOrder());
}

void
GlobalProblem::
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

        // hdlrAlgorithm->setInitialCondition(M_assembler->getInitialCondition());

        solveTimestep(t0, dt);

        M_assembler->setTimeAndPrevSolution(t0, hdlrAlgorithm->getSolution());
        M_assembler->postProcess();

        M_assembler->exportSolutions(t0, hdlrAlgorithm->getSolution());

        if (M_exportNorms)
        {
            M_assembler->appendNormsToFile(t0, hdlrAlgorithm->getSolution(),
                                                outFile);
        }
        if (M_exportErrors)
        {
            M_assembler->appendErrorsToFile(t0, hdlrAlgorithm->getSolution(),
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
        hdlrAlgorithm->setInitialCondition(M_assembler->getInitialCondition());
        M_assembler->exportSolutions(t, hdlrAlgorithm->getSolution());
        if (M_exportErrors)
        {
            M_assembler->appendErrorsToFile(t, hdlrAlgorithm->getSolution(),
                                                 outFile);
        }

        unsigned int count = 1;
        while (T - t > dt/2)
        {
            M_assembler->setTimestep(dt);
            solveTimestep(t, dt);
            t += dt;
            // we do this so that the single assemblers are aware of the
            // solution (e.g. for post processing)
            M_assembler->setTimeAndPrevSolution(t,
                                                     hdlrAlgorithm->getSolution());
            M_assembler->postProcess();
            if (count % save_every == 0)
                M_assembler->exportSolutions(t, hdlrAlgorithm->getSolution());
            if (M_exportNorms)
            {
                M_assembler->appendNormsToFile(t, hdlrAlgorithm->getSolution(),
                                                    outFile);
            }
            if (M_exportErrors)
            {
                M_assembler->appendErrorsToFile(t, hdlrAlgorithm->getSolution(),
                                                     outFile);
            }
            count++;
        }
    }
    if (M_exportNorms || M_exportErrors)
        outFile.close();
}

void
GlobalProblem::
setExportNorms(std::string filename)
{
    M_exportNorms = true;
    M_filename = filename;
}

void
GlobalProblem::
setExportErrors(std::string filename)
{
    M_exportErrors = true;
    M_filename = filename;
}

void
GlobalProblem::
setForcingFunction(FunctionType forcingFunction,
                   FunctionType forcingFunctionDt)
{
    M_assembler->setForcingFunction(forcingFunction, forcingFunctionDt);
}

void
GlobalProblem::
solveTimestep(const double& time, double& dt)
{
    M_timeMarchingAlgorithm->solveTimestep(time, dt);
}

void
GlobalProblem::
setLawInflow(std::function<double(double)> maxLaw)
{
    M_assembler->setLawInflow(maxLaw);
}

void
GlobalProblem::
setExactSolution(AbstractFunctor* exactSolution)
{
    M_assembler->setExactSolution(exactSolution);
}

void
GlobalProblem::
setLawDtInflow(std::function<double(double)> maxLawDt)
{
    M_assembler->setLawDtInflow(maxLawDt);
}

void
GlobalProblem::
printMeshSize(std::string filename)
{
    M_assembler->printMeshSize(filename);
}

}  // namespace RedMA
