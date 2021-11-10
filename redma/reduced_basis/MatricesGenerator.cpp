#include "MatricesGenerator.hpp"

namespace RedMA
{

MatricesGenerator::
MatricesGenerator(const DataContainer& data, EPETRACOMM comm) :
  M_data(data),
  M_comm(comm)
{
    if (M_comm->MyPID() != 0)
        throw new Exception("MatricesGenerator does not support more than one proc");
}

void
MatricesGenerator::
generate()
{
    createDefaultAssemblers();

    std::string outdir = "matrices";

    fs::create_directory(outdir);

    // dump norms
    for (auto& meshas : M_meshASPairMap)
    {
        fs::create_directory(outdir + "/" + meshas.first);

        unsigned int nComponents = meshas.second.first->getNumComponents();

        bool bcs = true;
        for (unsigned int i = 0; i < nComponents; i++)
        {
            std::string normstr = outdir + "/" + meshas.first + "/norm" +
                                  std::to_string(i);
            auto nnorm = spcast<aAssemblerFE>(meshas.second.first)->getNorm(i,bcs);
            spcast<SparseMatrix>(nnorm)->dump(normstr);
        }

        bcs = false;
        for (unsigned int i = 0; i < nComponents; i++)
        {
            std::string normstr = outdir + "/" + meshas.first + "/norm" +
                                  std::to_string(i) + "_nobcs";
            auto nnorm = spcast<aAssemblerFE>(meshas.second.first)->getNorm(i,bcs);
            convert<SparseMatrix>(nnorm)->dump(normstr);
        }
    }

    // dump matrices for supremizers
    unsigned int field2augment = M_data("rb/offline/basis/primal_supremizers/field2augment", 0);
    unsigned int limitingfield = M_data("rb/offline/basis/primal_supremizers/limitingfield", 1);

    for (auto& meshas : M_meshASPairMap)
    {
        auto constraintMatrix = spcast<aAssemblerFE>(meshas.second.first)->getConstraintMatrix();
        convert<SparseMatrix>(constraintMatrix)->dump(outdir + "/" + meshas.first + "/primalConstraint");
    }

    field2augment = M_data("rb/offline/basis/dual_supremizers/field2augment", 0);

    for (auto& meshas : M_meshASPairMap)
    {
        InterfaceAssembler interfaceAssembler(M_data);

        shp<BuildingBlock> buildingBlock = meshas.second.first->getTreeNode()->M_block;

        // create list of faces
        std::vector<GeometricFace> faces = buildingBlock->getOutlets();
	faces.push_back(buildingBlock->getInlet(0));
	std::cout << typeid(faces).name() << std::endl;
        for (auto face : faces)
        {
            shp<BlockMatrix> constraintMatrixBlock(new BlockMatrix(0,0));
            shp<BlockMatrix> constraintMatrixDummyBlock(new BlockMatrix(0,0));
            interfaceAssembler.buildCouplingMatrices(meshas.second.first,
                                                     face,
                                                     constraintMatrixBlock,
                                                     constraintMatrixDummyBlock);
            // we assume that the first block is the one to be coupled
            // (as in interface assembler)
            auto constraintMatrix = spcast<SparseMatrix>(constraintMatrixBlock->block(0,0));
            constraintMatrix->dump(outdir + "/" + meshas.first + "/dualConstraint" + std::to_string(face.M_flag));
        }
    }
}

void
MatricesGenerator::
createDefaultAssemblers()
{
    // using namespace boost::filesystem;

    std::string snapshotsdir = M_data("rb/offline/snapshots/directory", "snapshots");
    
    if (!fs::exists(snapshotsdir))
        throw new Exception("Snapshots directory has not been generated yet!");

    std::string paramdir = snapshotsdir + "/param";
    unsigned int i = 0;
    // we loop over the folders with the parameters
    while (fs::exists(paramdir + std::to_string(i)))
    {
	fs::directory_iterator end_it;
        for (fs::directory_iterator it(paramdir + std::to_string(i)); it != end_it; it++)
        { 	
            if (is_directory(it->status()))
            {
                std::string paramDir = it->path().string();
              
                unsigned int dashpos = paramDir.find_last_of("/");
                std::string nameMesh = paramDir.substr(dashpos + 1);
                
		if (M_meshASPairMap.find(nameMesh) == M_meshASPairMap.end())
                {
		    shp<TreeNode> defTreeNode = generateDefaultTreeNode(nameMesh);
                    shp<AssemblerType> defAssembler = AssemblerFactory(M_data, defTreeNode);
                    defAssembler->setup();

                    M_meshASPairMap[nameMesh].first = defAssembler;
                    M_meshASPairMap[nameMesh].second.resize(defAssembler->getNumComponents());
                    M_bases[nameMesh].reset(new RBBases(M_data, M_comm));
                    M_bases[nameMesh]->setNumberOfFields(defAssembler->getNumComponents());
                    unsigned int indexField = 0;
                    while (defAssembler->getFEspace(indexField))
                    {
                        M_bases[nameMesh]->setFESpace(defAssembler->getFEspace(indexField),
                                                      indexField);
                        indexField++;
                    }
                }
            }
        }

        i++;
    }
    printlog(MAGENTA, "done\n", M_data.getVerbose());
}

shp<TreeNode>
MatricesGenerator::
generateDefaultTreeNode(const std::string& nameMesh)
{	 	
    if (nameMesh.find("tube") != std::string::npos)
        return generateDefaultTube(nameMesh);
    else if (nameMesh.find("bifurcation_symmetric") != std::string::npos)
        return generateDefaultSymmetricBifurcation(nameMesh); 
    else if (nameMesh.find("aorta") != std::string::npos && nameMesh.find("bif1") == std::string::npos)
        return generateDefaultAortaBifurcation0(nameMesh);
    else if (nameMesh.find("aortabif1") != std::string::npos)
	return generateDefaultAortaBifurcation1(nameMesh);
    else
    	{
        throw new Exception("MatricesGenerator: this branch must still be implemented");
    	}	

    return nullptr;
}

shp<TreeNode>
MatricesGenerator::
generateDefaultTube(const std::string& nameMesh)
{
    unsigned int diameter = std::atoi(nameMesh.substr(5,6).c_str());
    unsigned int length = std::atoi(nameMesh.substr(7,8).c_str());
    std::string refinement = nameMesh.substr(9);

    shp<Tube> defaultTube(new Tube(M_comm, refinement, false, diameter, length));
    defaultTube->readMesh();
    defaultTube->setDiscretizationMethod("fem");
    // not very general
    defaultTube->setAssemblerType("navierstokes");

    shp<TreeNode> treeNode(new TreeNode(defaultTube, 1234));

    return treeNode;
}

shp<TreeNode>
MatricesGenerator::
generateDefaultSymmetricBifurcation(const std::string& nameMesh)
{
    unsigned int alpha = std::atoi(nameMesh.substr(13,15).c_str());
    std::string refinement = nameMesh.substr(17);

    shp<BifurcationSymmetric> defaultBifurcation(new BifurcationSymmetric(M_comm, refinement, false, alpha));
    defaultBifurcation->readMesh();
    defaultBifurcation->setDiscretizationMethod("fem");
    // not very general
    defaultBifurcation->setAssemblerType("navierstokes");

    shp<TreeNode> treeNode(new TreeNode(defaultBifurcation, 1234));

    return treeNode;
}

shp<TreeNode>
MatricesGenerator::
generateDefaultAortaBifurcation0(const std::string& nameMesh)
{
    std::cout << nameMesh.substr(10) << std::endl;
    std::string refinement = "normal";

    shp<AortaBifurcation0> defaultBifurcation(new AortaBifurcation0(M_comm, refinement));
    defaultBifurcation->readMesh();
    defaultBifurcation->setDiscretizationMethod("fem");
    
    defaultBifurcation->setAssemblerType("navierstokes");

    shp<TreeNode> treeNode(new TreeNode(defaultBifurcation, 1234));

    return treeNode;
}

shp<TreeNode>
MatricesGenerator::
generateDefaultAortaBifurcation1(const std::string& nameMesh)
{
    std::string refinement = "normal";
    std::cout << nameMesh.substr(10) << std::endl;
    shp<AortaBifurcation1> defaultBifurcation(new AortaBifurcation1(M_comm, refinement));
    defaultBifurcation->readMesh();
    defaultBifurcation->setDiscretizationMethod("fem");
    
    defaultBifurcation->setAssemblerType("navierstokes");

    shp<TreeNode> treeNode(new TreeNode(defaultBifurcation, 1234));

    return treeNode;
}

}// namespace RedMA
