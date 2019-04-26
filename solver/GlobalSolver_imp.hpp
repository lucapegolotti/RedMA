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

    std::string geometriesDir = datafile("geometric_structure/geometies_dir",
                                         "../../geometries/");

    M_timeMarchingAlgorithm =
            TimeMarchingAlgorithmsFactory<AssemblerType>(datafile);

    M_tree.readMeshes(geometriesDir);
    M_tree.traverseAndDeformGeometries();

    M_mapVector.reset(new MapVector());
    M_globalAssembler.buildPrimalStructures(M_tree, M_mapVector, M_globalMatrix);
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

    while (t < T)
    {
        solveTimestep(t);
        t += dt;
    }
}

template <class AssemblerType>
void
GlobalSolver<AssemblerType>::
solveTimestep(const double& time, double& dt)
{
    M_timeMarchingAlgorithm->solveTimestep(time, dt, M_globalAssembler);
}

}  // namespace RedMA
