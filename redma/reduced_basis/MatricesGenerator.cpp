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

    std::string assemblerType = M_data("assembler/type", "navierstokes");  // assuming the same for all

    // dump norm matrices (w\ and w\out BCs)
    for (auto& meshas : M_meshASPairMap)
    {
        fs::create_directory(outdir + "/" + meshas.first);

        unsigned int nComponents = meshas.second.first->getNumComponents();

        for (unsigned int i = 0; i < nComponents; i++)
        {
            std::string normstr_nobcs = outdir + "/" + meshas.first + "/norm" + std::to_string(i);
            auto nnorm_nobcs = spcast<aAssemblerFE>(meshas.second.first)->getNorm(i, false);
            spcast<SparseMatrix>(nnorm_nobcs)->dump(normstr_nobcs);

            if (i == meshas.second.first->getComponentBCs())
            {
                std::string normstr = outdir + "/" + meshas.first + "/norm" + std::to_string(i) + "_bcs";
                auto nnorm = spcast<aAssemblerFE>(meshas.second.first)->getNorm(i, true);
                spcast<SparseMatrix>(nnorm)->dump(normstr);
            }
        }

        if (!(std::strcmp(assemblerType.c_str(), "navierstokes_membrane")))
        {
            unsigned int i = nComponents;
            std::string normstr = outdir + "/" + meshas.first + "/norm" + std::to_string(i);
            auto nnorm = spcast<aAssemblerFE>(meshas.second.first)->getNorm(i, false);
            spcast<SparseMatrix>(nnorm)->dump(normstr);

        }
    }

    // dump velocity matrices
    for (auto& meshas : M_meshASPairMap)
    {
        auto massMatrix = spcast<aAssemblerFE>(meshas.second.first)->assembleMatrix(0);
        convert<SparseMatrix>(massMatrix->block(0,0))->dump(outdir + "/" + meshas.first + "/M");

        auto stiffnessMatrix = spcast<aAssemblerFE>(meshas.second.first)->assembleMatrix(1);
        convert<SparseMatrix>(stiffnessMatrix->block(0,0))->dump(outdir + "/" + meshas.first + "/A");
    }

    // dump pressure matrices
    for (auto& meshas : M_meshASPairMap)
    {
        auto divergenceMatrix = spcast<aAssemblerFE>(meshas.second.first)->assembleMatrix(2);
        convert<SparseMatrix>(divergenceMatrix->block(0,1))->dump(outdir + "/" + meshas.first + "/BdivT");
        convert<SparseMatrix>(divergenceMatrix->block(1,0))->dump(outdir + "/" + meshas.first + "/Bdiv");
    }

    // dump matrices for supremizers
    /*for (auto& meshas : M_meshASPairMap)
    {
        auto constraintMatrix = spcast<aAssemblerFE>(meshas.second.first)->getConstraintMatrix();
        convert<SparseMatrix>(constraintMatrix)->dump(outdir + "/" + meshas.first + "/primalConstraint");
    }*/

    for (auto& meshas : M_meshASPairMap)
    {
        InterfaceAssembler interfaceAssembler(M_data);
        InterfaceAssembler inletInflowAssembler(M_data);
        inletInflowAssembler.setAsInlet();
        InterfaceAssembler outletOutflowAssembler(M_data);
        outletOutflowAssembler.setAsOutlet();

        shp<BuildingBlock> buildingBlock = meshas.second.first->getTreeNode()->M_block;

        // create list of faces
        std::vector<GeometricFace> in_faces = buildingBlock->getInlets();
        std::vector<GeometricFace> out_faces = buildingBlock->getOutlets();
        std::vector<GeometricFace> faces;
	    faces.insert(std::end(faces), std::begin(in_faces), std::end(in_faces));
	    faces.insert(std::end(faces), std::begin(out_faces), std::end(out_faces));

	    // assemble coupling matrices
	    unsigned int cnt = 0;
        for (auto face : faces)
        {
            // Interface matrices
            shp<BlockMatrix> constraintMatrixBlockT(new BlockMatrix(0,0));
            shp<BlockMatrix> constraintMatrixBlock(new BlockMatrix(0,0));
            interfaceAssembler.buildCouplingMatrices(meshas.second.first,
                                                     face,
                                                     constraintMatrixBlockT,
                                                     constraintMatrixBlock);

            auto constraintMatrix = spcast<SparseMatrix>(constraintMatrixBlock->block(0,0));
            auto constraintMatrixT = spcast<SparseMatrix>(constraintMatrixBlockT->block(0,0));

            constraintMatrix->dump(outdir + "/" + meshas.first + "/B" + std::to_string(face.M_flag));
            constraintMatrixT->dump(outdir + "/" + meshas.first + "/BT" + std::to_string(face.M_flag));

            // constraintMatrix->dump(outdir + "/" + meshas.first + "/dualConstraint" + std::to_string(face.M_flag));

            // Weak inflow Dirichlet matrices
            constraintMatrixBlockT.reset(new BlockMatrix(0,0));
            constraintMatrixBlock.reset(new BlockMatrix(0,0));
            if (cnt < in_faces.size())
            {
                inletInflowAssembler.buildCouplingMatrices(meshas.second.first,
                                                            face,
                                                            constraintMatrixBlockT,
                                                            constraintMatrixBlock);

                constraintMatrix = spcast<SparseMatrix>(constraintMatrixBlock->block(0,0));
                constraintMatrixT = spcast<SparseMatrix>(constraintMatrixBlockT->block(0,0));

                constraintMatrix->dump(outdir + "/" + meshas.first + "/Bin" + std::to_string(face.M_flag));
                constraintMatrixT->dump(outdir + "/" + meshas.first + "/BinT" + std::to_string(face.M_flag));

                // constraintMatrix->dump(outdir + "/" + meshas.first + "/dualConstraintIn" + std::to_string(face.M_flag));
            }

            // Weak outflow Dirichlet matrices
            else
            {
                outletOutflowAssembler.buildCouplingMatrices(meshas.second.first,
                                                             face,
                                                             constraintMatrixBlockT,
                                                             constraintMatrixBlock);

                constraintMatrix = spcast<SparseMatrix>(constraintMatrixBlock->block(0,0));
                constraintMatrixT = spcast<SparseMatrix>(constraintMatrixBlockT->block(0,0));

                constraintMatrix->dump(outdir + "/" + meshas.first + "/Bout" + std::to_string(face.M_flag));
                constraintMatrixT->dump(outdir + "/" + meshas.first + "/BoutT" + std::to_string(face.M_flag));

                // constraintMatrix->dump(outdir + "/" + meshas.first + "/dualConstraintOut" + std::to_string(face.M_flag));
            }

            // TODO: export matrices for membrane assembler !!
            for (auto& meshas : M_meshASPairMap)
            {
                if (!(std::strcmp(assemblerType.c_str(), "navierstokes_membrane")))
                {
                    auto bcManager = meshas.second.first->getBCManager();

                    auto boundaryMassMatrixBlock = spcast<MembraneAssemblerFE>(meshas.second.first)->assembleBoundaryMass(bcManager);
                    spcast<SparseMatrix>(boundaryMassMatrixBlock->block(0, 0))->dump(outdir + "/" + meshas.first + "/M_bd");

                    auto boundaryStiffnessMatrixBlocks = spcast<MembraneAssemblerFE>(meshas.second.first)->assembleBoundaryStiffnessTerms(bcManager);
                    for (int i=0; i<boundaryStiffnessMatrixBlocks.size(); i++)
                        spcast<SparseMatrix>(boundaryStiffnessMatrixBlocks[i]->block(0, 0))->dump(outdir + "/" + meshas.first + "/A_bd_" + std::to_string(i));
                }
            }

            cnt++;
        }
    }
}

void
MatricesGenerator::
setDefaultParameterValues()
{
    M_data.setValueDouble("fluid/density", 1.06);
    M_data.setValueDouble("fluid/viscosity", 0.035);

    M_data.setValueDouble("structure/density", 1.2);
    M_data.setValueDouble("structure/thickness", 0.1);
    M_data.setValueInt("structure/constant_thickness", 1);  // enforce constant thickness in RB
    M_data.setValueDouble("structure/young", 4e6);
    M_data.setValueDouble("structure/poisson", 0.5);
}

void
MatricesGenerator::
createDefaultAssemblers()
{
    std::string snapshotsdir = M_data("rb/offline/snapshots/directory", "snapshots");

    if (!fs::exists(snapshotsdir))
        throw new Exception("Snapshots directory has not been generated yet!");

    setDefaultParameterValues();

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
    else if (nameMesh.find("bif_sym") != std::string::npos)
        return generateDefaultSymmetricBifurcation(nameMesh); 
    else if (nameMesh.find("aortabif0") != std::string::npos)
        return generateDefaultAortaBifurcation0(nameMesh);
    else if (nameMesh.find("aortabif1") != std::string::npos)
	    return generateDefaultAortaBifurcation1(nameMesh);
    else if (nameMesh.find("bypass") != std::string::npos)
        return generateDefaultBypass(nameMesh);
    else
        throw new Exception("[MatricesGenerator]: this branch must still be implemented");
}

shp<TreeNode>
MatricesGenerator::
generateDefaultTube(const std::string& nameMesh)
{
    unsigned int diameter = std::atoi(nameMesh.substr(5,6).c_str());
    unsigned int length = std::atoi(nameMesh.substr(7,8).c_str());
    std::string refinement = nameMesh.substr(9);

    std::string assemblerType = M_data("assembler/type", "navierstokes");

    shp<Tube> defaultTube(new Tube(M_comm, refinement, false, diameter, length));
    defaultTube->readMesh();

    defaultTube->setDatafile(M_data);
    defaultTube->setDiscretizationMethod("fem");
    defaultTube->setAssemblerType(assemblerType);

    shp<TreeNode> treeNode(new TreeNode(defaultTube, 1234));

    return treeNode;
}

shp<TreeNode>
MatricesGenerator::
generateDefaultSymmetricBifurcation(const std::string& nameMesh)
{
    unsigned int alpha = std::atoi(nameMesh.substr(13,15).c_str());
    std::string refinement = nameMesh.substr(17);

    std::string assemblerType = M_data("assembler/type", "navierstokes");

    shp<BifurcationSymmetric> defaultBifurcation(new BifurcationSymmetric(M_comm, refinement, false, alpha));
    defaultBifurcation->readMesh();

    defaultBifurcation->setDatafile(M_data);
    defaultBifurcation->setDiscretizationMethod("fem");
    defaultBifurcation->setAssemblerType(assemblerType);

    shp<TreeNode> treeNode(new TreeNode(defaultBifurcation, 1234));

    return treeNode;
}

shp<TreeNode>
MatricesGenerator::
generateDefaultAortaBifurcation0(const std::string& nameMesh)
{
    std::string refinement = "normal";

    std::string assemblerType = M_data("assembler/type", "navierstokes");

    shp<AortaBifurcation0> defaultBifurcation(new AortaBifurcation0(M_comm, refinement));
    defaultBifurcation->readMesh();

    defaultBifurcation->setDatafile(M_data);
    defaultBifurcation->setDiscretizationMethod("fem");
    defaultBifurcation->setAssemblerType(assemblerType);

    shp<TreeNode> treeNode(new TreeNode(defaultBifurcation, 1234));

    return treeNode;
}

shp<TreeNode>
MatricesGenerator::
generateDefaultAortaBifurcation1(const std::string& nameMesh)
{
    std::string refinement = "normal";

    std::string assemblerType = M_data("assembler/type", "navierstokes");

    shp<AortaBifurcation1> defaultBifurcation(new AortaBifurcation1(M_comm, refinement));
    defaultBifurcation->readMesh();

    defaultBifurcation->setDatafile(M_data);
    defaultBifurcation->setDiscretizationMethod("fem");
    defaultBifurcation->setAssemblerType(assemblerType);

    shp<TreeNode> treeNode(new TreeNode(defaultBifurcation, 1234));

    return treeNode;
}

shp<TreeNode>
MatricesGenerator::
generateDefaultBypass(const std::string &nameMesh)
{
    std::string assemblerType = M_data("assembler/type", "navierstokes");

    shp<Bypass> defaultBypass(new Bypass(M_comm));
    defaultBypass->readMesh();

    defaultBypass->setDatafile(M_data);
    defaultBypass->setDiscretizationMethod("fem");
    defaultBypass->setAssemblerType(assemblerType);

    shp<TreeNode> treeNode(new TreeNode(defaultBypass, 1234));

    return treeNode;
}

}// namespace RedMA
