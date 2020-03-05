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
    projectMDEIM();
    checkMDEIM();
    dumpMDEIMstructures();
}

void
MDEIMGenerator::
dumpMDEIMstructures()
{
    using namespace boost::filesystem;

    std::string outdir = M_data("rb/offline/mdeim/directory", "mdeims");

    create_directory(outdir);

    for (auto& mapit : M_primalBlockMDEIMsMap)
    {
        std::string curdir = outdir + "/" + mapit.first;
        create_directory(curdir);

        auto& mdeims = mapit.second;
        for (auto& mdeim : mdeims)
            mdeim.dumpMDEIMs(curdir);
    }
}

void
MDEIMGenerator::
checkMDEIM()
{
    for (auto& blockmdeims : M_primalBlockMDEIMsMap)
    {
        for (auto& mdeim : blockmdeims.second)
        {
            mdeim.prepareOnline();
            mdeim.checkOnSnapshots();
        }
    }

}

void
MDEIMGenerator::
takeMatricesSnapshots()
{
    using namespace boost::filesystem;

    unsigned int nSnapshots = M_data("rb/offline/mdeim/numbersnapshots", 50);
    double bound = M_data("rb/offline/mdeim/bound", 0.2);

    std::string outdir = M_data("rb/offline/mdeim/directory", "mdeims");

    if (exists(outdir))
        throw new Exception("Mdeims directory already exists!");

    for (unsigned int isnapshot = 0; isnapshot < nSnapshots; isnapshot++)
    {
        ProblemFEM problem(M_data, M_comm, false);
        problem.doStoreSolutions();

        problem.getTree().randomSampleAroundOriginalValue(bound);
        problem.setup();

        auto IDmeshTypeMap = problem.getBlockAssembler()->getIDMeshTypeMap();

        // get primal assemblers
        auto primalAssemblers = problem.getBlockAssembler()->getAssemblersMap();

        for (auto as : primalAssemblers)
        {
            unsigned int matCount = 0;
            std::vector<BlockMatrix<MatrixEp>> matrices = as.second->getMatrices();
            unsigned int nmatrices = matrices.size();

            if (M_primalBlockMDEIMsMap.find(IDmeshTypeMap[as.first]) == M_primalBlockMDEIMsMap.end())
                M_primalBlockMDEIMsMap[IDmeshTypeMap[as.first]].resize(nmatrices);

            for (unsigned int i = 0; i < nmatrices; i++)
            {
                auto& curmdeim = M_primalBlockMDEIMsMap[IDmeshTypeMap[as.first]][i];

                if (isnapshot == 0)
                {
                    curmdeim.setComm(M_comm);
                    curmdeim.setDataContainer(M_data);
                    curmdeim.resize(as.second->getNumComponents(),
                                    as.second->getNumComponents());
                    unsigned int fieldIndex = 0;
                    while (as.second->getFEspace(fieldIndex))
                    {
                        curmdeim.setFESpace(as.second->getFEspace(fieldIndex), fieldIndex);
                        curmdeim.setRangeMap(as.second->getFEspace(fieldIndex)->mapPtr(), fieldIndex);
                        curmdeim.setDomainMap(as.second->getFEspace(fieldIndex)->mapPtr(), fieldIndex);
                        fieldIndex++;
                    }
                    curmdeim.setMatrixIndex(i);
                    curmdeim.initialize(matrices[i]);
                }
                curmdeim.addSnapshot(matrices[i]);
            }
        }
    }
}

void
MDEIMGenerator::
performMDEIM()
{
    using namespace boost::filesystem;

    std::string outdir = M_data("rb/offline/mdeim/directory", "mdeims");
    create_directory(outdir);

    for (auto& mapit : M_primalBlockMDEIMsMap)
    {
        std::string curdir = outdir + "/" + mapit.first;
        create_directory(curdir);

        auto& mdeims = mapit.second;
        for (auto& mdeim : mdeims)
            mdeim.performMDEIM(curdir);
    }
}

void
MDEIMGenerator::
projectMDEIM()
{
    using namespace boost::filesystem;

    std::string basisdir = M_data("rb/offline/basis/directory", "basis");

    if (exists(basisdir))
    {
        for (auto& mdeims : M_primalBlockMDEIMsMap)
        {
            for (auto& mdeim : mdeims.second)
            {
                unsigned int numFields = mdeim.getNumRows();

                SHP(RBBases) curBases(new RBBases(M_data, M_comm));
                curBases->setPath(basisdir + "/" + mdeims.first);
                curBases->setNumberOfFields(numFields);

                for (unsigned int i = 0; i < numFields; i++)
                    curBases->setFESpace(mdeim.getFESpace(i), i);

                curBases->loadBases();
                mdeim.setRBBases(curBases);
                mdeim.projectMDEIMs();
            }
        }
    }
}

}  // namespace RedMA
