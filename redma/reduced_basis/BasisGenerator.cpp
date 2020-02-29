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

        std::string outdir = M_data("rbbasis/directory", "basis");

    if (boost::filesystem::exists(outdir))
        throw new Exception("Basis directory already exists!");

}

void
BasisGenerator::
generateBasis()
{
    parseFiles();
    performPOD();
    dumpBasis();
}

void
BasisGenerator::
parseFiles()
{
    printlog(MAGENTA, "[BasisGenerator] parsing files ... \n", M_data.getVerbose());

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
    printlog(MAGENTA, "done\n", M_data.getVerbose());
}

void
BasisGenerator::
performPOD()
{
    using namespace rbLifeV;

    double podtol = M_data("rbbasis/podtol", 1e-5);

    std::string outdir = M_data("rbbasis/directory", "basis");
    boost::filesystem::create_directory(outdir);

    printlog(MAGENTA, "[BasisGenerator] performing POD ... \n", M_data.getVerbose());
    for (auto pair : M_meshASPairMap)
    {
        unsigned int count = 0;
        VectorFunctions newBasisFunctions(pair.second.second.size());
        boost::filesystem::create_directory(outdir + "/" + pair.first);

        for (auto sn : pair.second.second)
        {
            ProperOrthogonalDecomposition pod(M_comm,
                                              pair.second.first->getFEspace(count)->map(),
                                              true);
            pod.initPOD(sn.size(), sn.data());
            pod.setSvdFileName(outdir + "/" + pair.first + "/svd" +
                               std::to_string(count) + ".txt");
            pod.generatePODbasisTol(podtol);

            unsigned int nbfs = pod.getRBdimension();
            std::vector<SHP(VECTOREPETRA)> basisFunctions(nbfs);

            pod.swapReducedBasis(basisFunctions, 0);
            newBasisFunctions[count] = basisFunctions;
            count++;
        }
        M_rbFunctions[pair.first] = newBasisFunctions;
    }
    printlog(MAGENTA, "done\n", M_data.getVerbose());
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
    infile.close();
}

void
BasisGenerator::
dumpBasis()
{
    using namespace boost::filesystem;

    std::string outdir = M_data("rbbasis/directory", "basis");
    bool binary = M_data("rbbasis/dumpbinary", true);

    std::ios_base::openmode omode = std::ios_base::app;
    if (binary)
        omode = omode | std::ios::binary;

    create_directory(outdir);

    for (auto meshbasis : M_rbFunctions)
    {
        std::string meshdir = outdir + "/" + meshbasis.first;
        create_directory(meshdir);

        unsigned int fieldIndex = 0;
        for (auto singlebasis : meshbasis.second)
        {
            std::ofstream outfile;
            outfile.open(meshdir + "/field" + std::to_string(fieldIndex) + ".basis", omode);

            for (auto vec : singlebasis)
            {
                FEVECTOR curvec;
                curvec.data() = vec;
                std::string str2write = curvec.getString(',') + "\n";

                outfile.write(str2write.c_str(), str2write.size());
            }

            outfile.close();
            fieldIndex++;
        }
    }

}

}  // namespace RedMA
