#include "MatricesGeneratorFixedGeometry.hpp"

namespace RedMA
{

MatricesGeneratorFixedGeometry::
MatricesGeneratorFixedGeometry(const DataContainer& data, EPETRACOMM comm) :
    aMatricesGenerator(data, comm)
{
}

void
MatricesGeneratorFixedGeometry::
generate()
{
    setDummyFlows();  // setting all inflows and outflows to 1
    createAssemblers();  // initializing the assembler

    std::string outdir = "matrices";
    fs::create_directory(outdir);

    std::string assemblerType = M_data("assembler/type", "navierstokes");

    unsigned int nComponents = M_assembler->getNumComponents();

    printlog(YELLOW, "[MatricesGeneratorFixedGeometry] Assembling norm matrices \n",
             this->M_data.getVerbose());

    for (unsigned int i = 0; i < nComponents; i++)
    {
        std::string normstr_nobcs = outdir + "/norm" + std::to_string(i);
        auto nnorm_nobcs = spcast<aAssemblerFE>(M_assembler)->getNorm(i, false);
        spcast<SparseMatrix>(nnorm_nobcs)->dump(normstr_nobcs);

        if (i == M_assembler->getComponentBCs())
        {
            std::string normstr = outdir + "/norm" + std::to_string(i) + "_bcs";
            auto nnorm = spcast<aAssemblerFE>(M_assembler)->getNorm(i, true);
            spcast<SparseMatrix>(nnorm)->dump(normstr);
        }
    }

    if (!(std::strcmp(assemblerType.c_str(), "navierstokes_membrane")))
    {
        unsigned int i = nComponents;
        std::string normstr = outdir + "/norm" + std::to_string(i);
        auto nnorm = spcast<aAssemblerFE>(M_assembler)->getNorm(i, false);
        spcast<SparseMatrix>(nnorm)->dump(normstr);

    }

    printlog(YELLOW, "[MatricesGeneratorFixedGeometry] Assembling primal matrices \n",
             this->M_data.getVerbose());

    auto bcManager = M_assembler->getBCManager();
    if (!(std::strcmp(assemblerType.c_str(), "navierstokes_membrane")))
        spcast<MembraneAssemblerFE>(M_assembler)->setAddBoundaryTerms(false);

    auto massMatrix = spcast<StokesAssemblerFE>(M_assembler)->assembleMass(bcManager);
    convert<SparseMatrix>(massMatrix->block(0,0))->dump(outdir + "/M");

    auto stiffnessMatrix = spcast<StokesAssemblerFE>(M_assembler)->assembleStiffness(bcManager);
    convert<SparseMatrix>(stiffnessMatrix->block(0,0))->dump(outdir + "/A");

    auto divergenceMatrix = spcast<StokesAssemblerFE>(M_assembler)->assembleDivergence(bcManager);
    convert<SparseMatrix>(divergenceMatrix->block(0,1))->dump(outdir + "/BdivT");
    convert<SparseMatrix>(divergenceMatrix->block(1,0))->dump(outdir + "/" + "/Bdiv");

    unsigned int n_cloths = M_data("cloth/n_cloths", 0);
    if (n_cloths > 0)
    {
        for (unsigned int i=0; i<n_cloths; i++)
        {
            auto clothMatrix = spcast<StokesAssemblerFE>(M_assembler)->assembleSingleBloodClothMatrix(bcManager, i);
            convert<SparseMatrix>(clothMatrix->block(0,0))->dump(outdir + "/Mcloth" + std::to_string(i));
        }
    }

    printlog(YELLOW, "[MatricesGeneratorFixedGeometry] Assembling dual matrices \n",
             this->M_data.getVerbose());

    Interface dummyInterfaceIn(nullptr, -1, M_assembler, 0, 0);
    Interface dummyInterfaceOut(M_assembler, 0, nullptr, -1, 0);

    shp<BuildingBlock> buildingBlock = M_assembler->getTreeNode()->M_block;
    bool hasNoSLipBCs = spcast<StokesAssemblerFE>(M_assembler)->hasNoSlipBCs();

    InterfaceAssembler interfaceAssembler(M_data, hasNoSLipBCs);
    InletInflowAssembler inletInflowAssembler(M_data, dummyInterfaceIn, hasNoSLipBCs, false);
    OutletOutflowAssembler outletOutflowAssembler(M_data, dummyInterfaceOut, hasNoSLipBCs, false);

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
        interfaceAssembler.buildCouplingMatrices(M_assembler,
                                                 face,
                                                 constraintMatrixBlockT,
                                                 constraintMatrixBlock);

        auto constraintMatrix = spcast<SparseMatrix>(constraintMatrixBlock->block(0,0));
        auto constraintMatrixT = spcast<SparseMatrix>(constraintMatrixBlockT->block(0,0));

        constraintMatrix->dump(outdir + "/B" + std::to_string(face.M_flag));
        constraintMatrixT->dump(outdir + "/BT" + std::to_string(face.M_flag));

        // Weak inflow Dirichlet matrices
        constraintMatrixBlockT.reset(new BlockMatrix(0,0));
        constraintMatrixBlock.reset(new BlockMatrix(0,0));
        if (cnt < in_faces.size())
        {
            inletInflowAssembler.buildCouplingMatrices(M_assembler,
                                                       face,
                                                       constraintMatrixBlockT,
                                                       constraintMatrixBlock);

            constraintMatrix = spcast<SparseMatrix>(constraintMatrixBlock->block(0,0));
            constraintMatrixT = spcast<SparseMatrix>(constraintMatrixBlockT->block(0,0));

            constraintMatrix->multiplyByScalar(-1.0);
            constraintMatrixT->multiplyByScalar(-1.0);

            constraintMatrix->dump(outdir + "/Bin" + std::to_string(face.M_flag));
            constraintMatrixT->dump(outdir + "/BinT" + std::to_string(face.M_flag));

            shp<BlockVector> rhsVector(new BlockVector(0));
            inletInflowAssembler.buildRhsVector(M_assembler, face, rhsVector, 0.0);
            rhsVector->block(0)->multiplyByScalar(-1.0);
            rhsVector->block(0)->dump(outdir + "/RHS_in" + std::to_string(face.M_flag));
        }

        // Weak outflow Dirichlet matrices
        else
        {
            outletOutflowAssembler.buildCouplingMatrices(M_assembler,
                                                         face,
                                                         constraintMatrixBlockT,
                                                         constraintMatrixBlock);

            constraintMatrix = spcast<SparseMatrix>(constraintMatrixBlock->block(0,0));
            constraintMatrixT = spcast<SparseMatrix>(constraintMatrixBlockT->block(0,0));

            constraintMatrix->dump(outdir + "/Bout" + std::to_string(face.M_flag));
            constraintMatrixT->dump(outdir + "/BoutT" + std::to_string(face.M_flag));

            shp<BlockVector> rhsVector(new BlockVector(0));
            outletOutflowAssembler.buildRhsVector(M_assembler, face, rhsVector, 0.0);
            rhsVector->block(0)->dump(outdir + "/RHS_out" + std::to_string(face.M_flag));
        }

        // flow rate vectors at I/O
        shp<VECTOREPETRA> flowRateVector;
        flowRateVector.reset(new VECTOREPETRA(*(spcast<StokesAssemblerFE>(M_assembler)->getFlowRateVector(face.M_flag)),
                                              LifeV::Unique));
        std::string filename;
        if (cnt < in_faces.size())
            filename = outdir + "/q_in" + std::to_string(cnt);
        else
            filename = outdir + "/q_out" + std::to_string(cnt - in_faces.size());
        flowRateVector->spy(filename);

        // boundary matrices, if the membrane model is selected
        if (!(std::strcmp(assemblerType.c_str(), "navierstokes_membrane")))
        {
            auto bcManager = M_assembler->getBCManager();

            auto boundaryMassMatrixBlock = spcast<MembraneAssemblerFE>(M_assembler)->assembleBoundaryMass(bcManager);
            spcast<SparseMatrix>(boundaryMassMatrixBlock->block(0, 0))->dump(outdir + "/M_bd");

            auto boundaryStiffnessMatrixBlocks = spcast<MembraneAssemblerFE>(M_assembler)->assembleBoundaryStiffnessTerms(bcManager);
            for (int i=0; i<boundaryStiffnessMatrixBlocks.size(); i++)
                spcast<SparseMatrix>(boundaryStiffnessMatrixBlocks[i]->block(0, 0))->dump(outdir + "/A_bd_" + std::to_string(i));
        }

        cnt++;
    }
}

void
MatricesGeneratorFixedGeometry::
createAssemblers()
{
    std::string param_type = M_data("rb/offline/snapshots/param_type", "geometric");
    std::list<std::string> param_types = M_data.stringTokenizer(param_type, ',');

    if (std::find(std::begin(param_types), std::end(param_types), "physics") != std::end(param_types))
    {
        std::map<std::string, bool> categories;

        categories["fluid"] = M_data("rb/offline/snapshots/sample_fluid_physics", false);
        categories["structure"] = M_data("rb/offline/snapshots/sample_structure_physics", false);
        categories["cloth"] = M_data("rb/offline/snapshots/sample_cloth_physics", false);

        setDefaultParameterValues(categories);
    }

    shp<TreeNode> treeNode = generateTreeNode();

    M_assembler = AssemblerFactory(M_data, treeNode);
    M_assembler->setup();
}

shp<TreeNode>
MatricesGeneratorFixedGeometry::
generateTreeNode()
{
    GeometryParser geometryParser = GeometryParser(M_data,
                                                   M_data("geometric_structure/xmlfile","tree.xml"),
                                                   M_comm, M_data.getVerbose());
    TreeStructure tree = geometryParser.getTree();

    std::string geometriesDir = M_data("geometric_structure/geometries_dir",
                                       "../../../meshes/");
    tree.readMeshes(geometriesDir);
    tree.traverseAndDeformGeometries(true);

    shp<TreeNode> treeNode = tree.getRoot();

    return treeNode;
}

}// namespace RedMA

