#include "MDEIMGenerator.hpp"

namespace RedMA
{

MDEIMGenerator::
MDEIMGenerator(const DataContainer& data, EPETRACOMM comm) :
  M_data(data),
  M_comm(comm)
{
    if (M_comm->MyPID() != 0)
        throw new Exception("MDEIMGenerator does not support more than one proc");
}

void
MDEIMGenerator::
generateMDEIM()
{
    takeMatricesSnapshots();
    performMDEIM();
    checkMDEIM();
}

void
MDEIMGenerator::
checkMDEIM()
{
    unsigned int nChecks = M_data("mdeim/checksonline", 5);
    double bound = M_data("mdeim/bound", 0.2);

    for (unsigned int i = 0; i < nChecks; i++)
    {
        ProblemFEM problem(M_data, M_comm, false);

        problem.getTree().randomSampleAroundOriginalValue(bound);
        problem.setup();

        auto assemblers = problem.getBlockAssembler()->getAssemblersMap();
        auto IDmeshTypeMap = problem.getBlockAssembler()->getIDMeshTypeMap();

        for (auto as : assemblers)
        {
            for (auto& mdeim : M_mdeimsMap[IDmeshTypeMap[as.first]])
            {
                mdeim.setAssembler(as.second);
                mdeim.prepareOnline();
                mdeim.checkOnline();
                // mdeim.checkOnSnapshots();
            }
        }
    }

}

void
MDEIMGenerator::
takeMatricesSnapshots()
{
    unsigned int nSnapshots = M_data("mdeim/numbersnapshots", 50);
    double bound = M_data("mdeim/bound", 0.2);
    int seed = M_data("mdeim/seed", 1234);

    for (unsigned int i = 0; i < nSnapshots; i++)
    {
        ProblemFEM problem(M_data, M_comm, false);

        problem.getTree().randomSampleAroundOriginalValue(bound);
        problem.setup();

        auto assemblers = problem.getBlockAssembler()->getAssemblersMap();
        auto IDmeshTypeMap = problem.getBlockAssembler()->getIDMeshTypeMap();

        for (auto as : assemblers)
        {
            unsigned int matCount = 0;
            std::vector<BlockMatrix<MatrixEp>> matrices = as.second->getMatrices();
            unsigned int nmatrices = matrices.size();

            if (M_mdeimsMap.find(IDmeshTypeMap[as.first]) == M_mdeimsMap.end())
                M_mdeimsMap[IDmeshTypeMap[as.first]].resize(nmatrices);

            for (unsigned int i = 0; i < nmatrices; i++)
            {
                M_mdeimsMap[IDmeshTypeMap[as.first]][i].setComm(M_comm);
                M_mdeimsMap[IDmeshTypeMap[as.first]][i].setDataContainer(M_data);
                M_mdeimsMap[IDmeshTypeMap[as.first]][i].setAssembler(assemblers[as.first]);
                M_mdeimsMap[IDmeshTypeMap[as.first]][i].addSnapshot(matrices[i]);
                M_mdeimsMap[IDmeshTypeMap[as.first]][i].setMatrixIndex(i);
            }
        }
    }
}

void
MDEIMGenerator::
performMDEIM()
{
    for (auto& mapit : M_mdeimsMap)
    {
        auto& mdeims = mapit.second;
        for (auto& mdeim : mdeims)
            mdeim.performMDEIM();
    }
}



}  // namespace RedMA
