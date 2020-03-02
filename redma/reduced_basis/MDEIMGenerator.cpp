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

    std::string outdir = M_data("mdeim/directory", "mdeims");

    create_directory(outdir);

    for (auto& mapit : M_blockMDEIMsMap)
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
    for (auto& blockmdeims : M_blockMDEIMsMap)
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

    unsigned int nSnapshots = M_data("mdeim/numbersnapshots", 50);
    double bound = M_data("mdeim/bound", 0.2);

    std::string outdir = M_data("mdeim/directory", "mdeims");

    if (exists(outdir))
        throw new Exception("Mdeims directory already exists!");

    for (unsigned int i = 0; i < nSnapshots; i++)
    {
        ProblemFEM problem(M_data, M_comm, false);
        problem.doStoreSolutions();

        problem.getTree().randomSampleAroundOriginalValue(bound);
        problem.setup();

        auto assemblers = problem.getBlockAssembler()->getAssemblersMap();
        auto IDmeshTypeMap = problem.getBlockAssembler()->getIDMeshTypeMap();

        for (auto as : assemblers)
        {
            unsigned int matCount = 0;
            std::vector<BlockMatrix<MatrixEp>> matrices = as.second->getMatrices();
            unsigned int nmatrices = matrices.size();

            if (M_blockMDEIMsMap.find(IDmeshTypeMap[as.first]) == M_blockMDEIMsMap.end())
                M_blockMDEIMsMap[IDmeshTypeMap[as.first]].resize(nmatrices);

            for (unsigned int i = 0; i < nmatrices; i++)
            {
                M_blockMDEIMsMap[IDmeshTypeMap[as.first]][i].setComm(M_comm);
                M_blockMDEIMsMap[IDmeshTypeMap[as.first]][i].setDataContainer(M_data);
                M_blockMDEIMsMap[IDmeshTypeMap[as.first]][i].setAssembler(assemblers[as.first]);
                M_blockMDEIMsMap[IDmeshTypeMap[as.first]][i].addSnapshot(matrices[i]);
                M_blockMDEIMsMap[IDmeshTypeMap[as.first]][i].setMatrixIndex(i);
            }
        }
    }
}

void
MDEIMGenerator::
performMDEIM()
{
    using namespace boost::filesystem;

    std::string outdir = M_data("mdeim/directory", "mdeims");
    create_directory(outdir);

    for (auto& mapit : M_blockMDEIMsMap)
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

    std::string basisdir = M_data("rbbasis/directory", "basis");

    if (exists(basisdir))
    {
        for (auto& mdeims : M_blockMDEIMsMap)
        {
            for (auto& mdeim : mdeims.second)
            {
                unsigned int numFields = mdeim.getAssembler()->getNumComponents();

                SHP(RBBases) curBases(new RBBases(M_data, M_comm));
                curBases->setPath(basisdir + "/" + mdeims.first);
                curBases->setNumberOfFields(numFields);

                for (unsigned int i = 0; i < numFields; i++)
                    curBases->setFESpace(mdeim.getAssembler()->getFEspace(i), i);

                curBases->loadBases();
                mdeim.setRBBases(curBases);
                mdeim.projectMDEIMs();
            }
        }
    }
}

std::vector<SHP(VECTOREPETRA)>
MDEIMGenerator::
readRBBasisFromFile(std::string file, SHP(FESPACE) fespace)
{
    std::vector<SHP(VECTOREPETRA)> retBasis;

    std::ifstream infile(file);
    std::string line;
    while(std::getline(infile,line))
    {
        SHP(VECTOREPETRA) newVector(new VECTOREPETRA(fespace->map()));

        std::stringstream linestream(line);
        std::string value;
        unsigned int i = 0;
        while(getline(linestream,value,','))
        {
            newVector->operator[](i) = std::atof(value.c_str());
            i++;
        }
        if (i != newVector->epetraVector().GlobalLength())
            throw new Exception("Stored basis length does not match fespace dimension!");

        retBasis.push_back(newVector);
    }
    infile.close();
    return retBasis;
}

}  // namespace RedMA
