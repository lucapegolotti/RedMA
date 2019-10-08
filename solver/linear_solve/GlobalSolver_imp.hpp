// implementation of template class

namespace RedMA
{

template <class AssemblerType>
GlobalSolver<AssemblerType>::
GlobalSolver(const GetPot& datafile, commPtr_Type comm, bool verbose) :
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

    M_globalAssembler.setup(M_tree);

    M_timeMarchingAlgorithm =
            TimeMarchingAlgorithmsFactory<AssemblerType>(datafile,
                                                         &M_globalAssembler,
                                                         M_comm, M_verbose);

    M_globalAssembler.setTimeIntegrationOrder(M_timeMarchingAlgorithm->getOrder());
}

template <class AssemblerType>
void
GlobalSolver<AssemblerType>::
solve()
{
    double t0 = M_datafile("time_discretization/t0", 0.0);
    double T = M_datafile("time_discretization/T", 1.0);
    double dt = M_datafile("time_discretization/dt", 0.01);
    unsigned int save_every = M_datafile("exporter/save_every", 1);

    double t = t0;
    TimeMarchingAlgorithmPtr hdlrAlgorithm = M_timeMarchingAlgorithm;
    hdlrAlgorithm->setInitialCondition(M_globalAssembler.getInitialCondition());
    M_globalAssembler.exportSolutions(t, hdlrAlgorithm->getSolution());
    std::ofstream outFile;
    if (M_exportNorms)
    {
        outFile.open(M_filename);
        outFile << AssemblerType::normFileFirstLine() << std::flush;
    }
    if (M_exportErrors)
    {
        outFile.open(M_filename);
        outFile << AssemblerType::errorsFileFirstLine() << std::flush;
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
    if (M_exportNorms || M_exportErrors)
        outFile.close();
}

template <class AssemblerType>
void
GlobalSolver<AssemblerType>::
setExportNorms(std::string filename)
{
    M_exportNorms = true;
    M_filename = filename;
}

template <class AssemblerType>
void
GlobalSolver<AssemblerType>::
setExportErrors(std::string filename)
{
    M_exportErrors = true;
    M_filename = filename;
}

template <class AssemblerType>
void
GlobalSolver<AssemblerType>::
solveTimestep(const double& time, double& dt)
{
    M_timeMarchingAlgorithm->solveTimestep(time, dt);
}

template <class AssemblerType>
void
GlobalSolver<AssemblerType>::
setLawInflow(std::function<double(double)> maxLaw)
{
    M_globalAssembler.setLawInflow(maxLaw);
}

template <class AssemblerType>
void
GlobalSolver<AssemblerType>::
setExactSolution(AbstractFunctor* exactSolution)
{
    M_globalAssembler.setExactSolution(exactSolution);
}

template <class AssemblerType>
void
GlobalSolver<AssemblerType>::
setLawDtInflow(std::function<double(double)> maxLawDt)
{
    M_globalAssembler.setLawDtInflow(maxLawDt);
}

}  // namespace RedMA
