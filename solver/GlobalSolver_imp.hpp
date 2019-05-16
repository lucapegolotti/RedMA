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
}

template <class AssemblerType>
void
GlobalSolver<AssemblerType>::
solve()
{
    double t0 = M_datafile("time_discretization/t0", 0.0);
    double T = M_datafile("time_discretization/T", 1.0);
    double dt = M_datafile("time_discretization/dt", 0.01);

    double t = t0;
    TimeMarchingAlgorithmPtr hdlrAlgorithm = M_timeMarchingAlgorithm;
    M_globalAssembler.exportSolutions(t, hdlrAlgorithm->getSolution());
    while (t < T)
    {
        solveTimestep(t, dt);
        t += dt;
        M_globalAssembler.exportSolutions(t, hdlrAlgorithm->getSolution());
    }
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
setMaxVelocityLawInflow(std::function<double(double)> maxLaw)
{
    M_globalAssembler.setMaxVelocityLawInflow(maxLaw);
}

template <class AssemblerType>
void
GlobalSolver<AssemblerType>::
setMaxVelocityDtLawInflow(std::function<double(double)> maxLawDt)
{
    M_globalAssembler.setMaxVelocityDtLawInflow(maxLawDt);
}

}  // namespace RedMA
