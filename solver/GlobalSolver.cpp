#include <GlobalSolver.hpp>

namespace RedMA
{

template<class AssemblerType>
GlobalSolver<AssemblerType>::
GlobalSolver(std::string xmlFile, std::string geometriesDir,
             commPtr_Type comm, bool verbose) :
  M_geometryParser(xmlFile, comm, verbose)
{
    M_tree = M_geometryParser.getTree();

    M_tree.readMeshes(geometriesDir);
    M_tree.traverseAndDeformGeometries();

    GlobalAssembler<AssemblerType> globalAssembler;

    globalAssembler.buildPrimalStructures(M_tree, M_mapVector, M_globalMatrix);
}

}  // namespace RedMA
