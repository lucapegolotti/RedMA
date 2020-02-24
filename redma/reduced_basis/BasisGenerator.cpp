#include "BasisGenerator.hpp"

namespace RedMA
{

BasisGenerator::
BasisGenerator(const DataContainer& data, EPETRACOMM comm) :
  M_data(data),
  M_comm(comm)
{
    if (M_comm->MyPID() != 0)
        throw new Exception("BasisGenerator does not support more than one proc");
}

void
BasisGenerator::
generateBasis()
{
    using namespace boost::filesystem;

    std::string snapshotsdir = M_data("snapshots/directory", "snapshots");

    if (!exists(snapshotsdir))
        throw new Exception("Snapshots directory has not been generated yet!");

    std::string paramdir = snapshotsdir + "/param";
    unsigned int i = 0;
    // we loop over the folders with the parameters
    while (exists(paramdir + std::to_string(i)))
    {
        directory_iterator end_it;
        for (directory_iterator it(paramdir + std::to_string(i)); it != end_it; it++)
        {
            if (is_directory(it->status()))
                parseParameterSnapshots(it->path().string());
        }

        i++;
    }
}

SHP(TreeNode)
BasisGenerator::
generateDefaultTreeNode(const std::string& nameMesh)
{
    if (nameMesh.find("tube") != std::string::npos)
        return generateDefaultTube(nameMesh);
    else if (nameMesh.find("bifurcation_symmetric"))
    {
        throw new Exception("BasisGenerator: this branch must still be implemented");
    }

    return nullptr;
}

SHP(TreeNode)
BasisGenerator::
generateDefaultTube(const std::string& nameMesh)
{
    unsigned int diameter = std::atoi(nameMesh.substr(5,6).c_str());
    unsigned int length = std::atoi(nameMesh.substr(7,8).c_str());
    std::string refinement = nameMesh.substr(9);

    SHP(Tube) defaultTube(new Tube(M_comm, refinement, false, diameter, length));
    defaultTube->readMesh();

    SHP(TreeNode) treeNode(new TreeNode(defaultTube, 1234));

    return treeNode;
}

void
BasisGenerator::
parseParameterSnapshots(const std::string& paramDir)
{
    using namespace boost::filesystem;

    unsigned int dashpos = paramDir.find_last_of("/");
    std::string nameMesh = paramDir.substr(dashpos + 1);

    if (M_meshASPairMap.find(nameMesh) == M_meshASPairMap.end())
    {
        SHP(TreeNode) defTreeNode = generateDefaultTreeNode(nameMesh);
        SHP(AssemblerType) defAssembler = AssemblerFactory<FEVECTOR COMMA FEMATRIX>(M_data, defTreeNode);
        defAssembler->setup();

        M_meshASPairMap[nameMesh].first = defAssembler;
        M_meshASPairMap[nameMesh].second.resize(defAssembler->getNumComponents());
    }

    directory_iterator end_it;
    for (directory_iterator it(paramDir); it != end_it; it++)
    {
        std::string curfile = it->path().string();
        if (curfile.find(".snap") != std::string::npos)
        {
            unsigned int dotpos = curfile.find_last_of(".");
            unsigned int componentIndex = std::atoi(curfile.substr(dotpos-1,dotpos).c_str());

            addSnapshotsFromFile(curfile,
                                 M_meshASPairMap[nameMesh].second[componentIndex],
                                 M_meshASPairMap[nameMesh].first->getFEspace(componentIndex));
        }
    }
}

void
BasisGenerator::
addSnapshotsFromFile(const std::string& snapshotsFile,
                     std::vector<SHP(VECTOREPETRA)>& snapshots,
                     SHP(FESPACE) fespace)
{

    std::ifstream infile(snapshotsFile);
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
            throw new Exception("Stored snapshot length does not match fespace dimension!");

        snapshots.push_back(newVector);
    }
}

}  // namespace RedMA
