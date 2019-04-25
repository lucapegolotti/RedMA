// implementation of template class

namespace RedMA
{

template <class AssemblerType>
GlobalSolver<AssemblerType>::
GlobalSolver(GetPot datafile, commPtr_Type comm, bool verbose) :
  M_geometryParser(datafile("geometric_structure/xmlfile","tree.xml"),
                   comm, verbose),
  M_datafile(datafile)
{
    M_tree = M_geometryParser.getTree();

    std::string geometriesDir = datafile("geometric_structure/geometies_dir",
                                         "../../geometries/");

    M_tree.readMeshes(geometriesDir);
    M_tree.traverseAndDeformGeometries();

    GlobalAssembler<AssemblerType> globalAssembler;

    globalAssembler.buildPrimalStructures(M_tree, M_mapVector, M_globalMatrix);
}

}  // namespace RedMA
