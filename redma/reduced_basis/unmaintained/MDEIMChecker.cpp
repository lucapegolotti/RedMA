#include "MDEIMChecker.hpp"

namespace RedMA
{

MDEIMChecker::
MDEIMChecker(const DataContainer& data, EPETRACOMM comm) :
  M_data(data),
  M_comm(comm)
{
    // if (M_comm->MyPID() != 0)
    //     throw new Exception("MDEIMChecker does not support more than one proc");
    //
    // if (!M_data("rb/online/mdeim/loadfullbasis", false))
    //     throw new Exception("rb/online/mdeim/loadfullbasis == false!");
    //
    // loadMDEIMs();
}

void
MDEIMChecker::
checkOnlineMDEIM()
{
    // unsigned int nChecks = M_data("rb/offline/mdeim/checksonline", 5);
    // double bound = M_data("rb/offline/mdeim/bound", 0.2);
    //
    // for (unsigned int i = 0; i < nChecks; i++)
    // {
    //     GlobalProblem problem(M_data, M_comm, false);
    //
    //     problem.getTree().randomSampleAroundOriginalValue(bound);
    //     // this can be optimized. Matrices are assembled twice
    //     problem.setup();
    //
    //     auto assemblers = problem.getBlockAssembler()->getAssemblersMap();
    //     auto IDmeshTypeMap = problem.getBlockAssembler()->getIDMeshTypeMap();
    //
    //     for (auto as : assemblers)
    //     {
    //         for (auto& mdeim : M_blockMDEIMsMap[IDmeshTypeMap[as.first]])
    //         {
    //             BlockMatrix<MatrixEp> fullMatrix = as.second->assembleMatrix(mdeim.getIndex());
    //             BlockMatrix<MatrixEp> reducedMatrix = as.second->assembleMatrix(mdeim.getIndex(),
    //                                                                             &mdeim.getMDEIMStructure());
    //             mdeim.checkOnline(fullMatrix, reducedMatrix);
    //         }
    //     }
    // }
}

void
MDEIMChecker::
loadMDEIMs()
{
    // // using namespace boost::filesystem;
    //
    // // this is only to get the assembler to give to the mdeims
    // GlobalProblem problem(M_data, M_comm, false);
    // problem.setup();
    // auto assemblers = problem.getBlockAssembler()->getAssemblersMap();
    // auto IDmeshTypeMap = problem.getBlockAssembler()->getIDMeshTypeMap();
    //
    // std::string outdir = M_data("rb/offline/mdeim/directory", "mdeims");
    //
    // if (!exists(outdir))
    //     throw new Exception("MDEIMS directory does not exist");
    //
    // directory_iterator end_it;
    // for (directory_iterator it(outdir); it != end_it; it++)
    // {
    //     std::string curdir = it->path().string();
    //     unsigned int dashpos = curdir.find_last_of("/");
    //     std::string meshName = curdir.substr(dashpos+1);
    //
    //     for (directory_iterator itIn(curdir); itIn != end_it; itIn++)
    //     {
    //         curdir = itIn->path().string();
    //
    //         int matrixIndex = std::atoi(curdir.substr(curdir.size()-1).c_str());
    //
    //         BlockMDEIM newMdeim;
    //         newMdeim.setComm(M_comm);
    //         newMdeim.setDataContainer(M_data);
    //
    //         unsigned int id;
    //         // wow, so ugly
    //         for (auto idmesh : IDmeshTypeMap)
    //             if (idmesh.second == meshName)
    //                 id = idmesh.first;
    //
    //         newMdeim.resize(assemblers[id]->getNumComponents(),
    //                         assemblers[id]->getNumComponents());
    //
    //         unsigned int fieldIndex = 0;
    //         while (assemblers[id]->getFEspace(fieldIndex))
    //         {
    //             newMdeim.setFESpace(assemblers[id]->getFEspace(fieldIndex), fieldIndex);
    //             newMdeim.setRangeMap(assemblers[id]->getFEspace(fieldIndex)->mapPtr(), fieldIndex);
    //             newMdeim.setDomainMap(assemblers[id]->getFEspace(fieldIndex)->mapPtr(), fieldIndex);
    //             fieldIndex++;
    //         }
    //
    //         newMdeim.setMatrixIndex(matrixIndex);
    //         newMdeim.loadMDEIM(curdir);
    //         M_blockMDEIMsMap[meshName].push_back(newMdeim);
    //     }
    // }

}

}  // namespace RedMA
