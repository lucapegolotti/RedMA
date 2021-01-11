#include "BasisGenerator.hpp"

namespace RedMA
{

BasisGenerator::
BasisGenerator(const DataContainer& data, EPETRACOMM comm) :
  M_data(data),
  M_comm(comm)
{
    // if (M_comm->MyPID() != 0)
    //     throw new Exception("BasisGenerator does not support more than one proc");
    //
    // std::string outdir = M_data("rb/offline/basis/directory", "basis");
    //
    // if (fs::esists(outdir))
    //     printlog(YELLOW,"[BasisGenerator] basis directory already exists!\n", data.getVerbose());
    //
    // // we want to consider the whole basis when adding supremizers
    // M_data.setValueDouble("rb/online/basis/podtol", 0.0);
}

void
BasisGenerator::
generateBasis()
{
    // createDefaultAssemblers();
    // parseFiles();
    // performPOD();
    // addSupremizers();
    // orthonormalize();
    // dumpBasis();
}

void
BasisGenerator::
generateMatricesOnly()
{
    // createDefaultAssemblers();
    //
    // std::string outdir = "matricesForOffline";
    //
    // fs::create_directory(outdir);
    //
    // // dump norms
    // for (auto& meshas : M_meshASPairMap)
    // {
    //     fs::create_directory(outdir + "/" + meshas.first);
    //
    //     unsigned int nComponents = meshas.second.first->getNumComponents();
    //
    //     bool bcs = true;
    //     for (unsigned int i = 0; i < nComponents; i++)
    //         meshas.second.first->getNorm(i,bcs).dump(outdir + "/" + meshas.first + "/norm" +
    //                              std::to_string(i));
    //
    //     bcs = false;
    //     for (unsigned int i = 0; i < nComponents; i++)
    //         meshas.second.first->getNorm(i,bcs).dump(outdir + "/" + meshas.first + "/norm" +
    //                              std::to_string(i) + "_nobcs");
    // }
    //
    // // dump matrices for supremizers
    // unsigned int field2augment = M_data("rb/offline/basis/primal_supremizers/field2augment", 0);
    // unsigned int limitingfield = M_data("rb/offline/basis/primal_supremizers/limitingfield", 1);
    //
    // for (auto& meshas : M_meshASPairMap)
    // {
    //     MatrixEp constraintMatrix = meshas.second.first->getConstraintMatrix();
    //     constraintMatrix.dump(outdir + "/" + meshas.first + "/primalConstraint");
    // }
    //
    // field2augment = M_data("rb/offline/basis/dual_supremizers/field2augment", 0);
    //
    // for (auto& meshas : M_meshASPairMap)
    // {
    //     InterfaceAssembler<VectorEp, MatrixEp> interfaceAssembler(M_data);
    //
    //     SHP(BuildingBlock) buildingBlock = meshas.second.first->getTreeNode()->M_block;
    //
    //     // create list of faces
    //     std::vector<GeometricFace> faces = buildingBlock->getOutlets();
    //     faces.push_back(buildingBlock->getInlet());
    //
    //     for (auto face : faces)
    //     {
    //         BlockMatrix<MatrixEp> constraintMatrixBlock;
    //         BlockMatrix<MatrixEp> constraintMatrixDummyBlock;
    //         interfaceAssembler.buildCouplingMatrices(meshas.second.first,
    //                                                  face,
    //                                                  constraintMatrixBlock,
    //                                                  constraintMatrixDummyBlock);
    //         // we assume that the first block is the one to be coupled
    //         // (as in interface assembler)
    //         MatrixEp constraintMatrix = constraintMatrixBlock.block(0,0);
    //         constraintMatrix.dump(outdir + "/" + meshas.first + "/dualConstraint" + std::to_string(face.M_flag));
    //     }
    // }
}

void
BasisGenerator::
createDefaultAssemblers()
{
    // // using namespace boost::filesystem;
    //
    // std::string snapshotsdir = M_data("rb/offline/snapshots/directory", "snapshots");
    //
    // if (!exists(snapshotsdir))
    //     throw new Exception("Snapshots directory has not been generated yet!");
    //
    // std::string paramdir = snapshotsdir + "/param";
    // unsigned int i = 0;
    // // we loop over the folders with the parameters
    // while (fs::exists(paramdir + std::to_string(i)))
    // {
    //     directory_iterator end_it;
    //     for (directory_iterator it(paramdir + std::to_string(i)); it != end_it; it++)
    //     {
    //         if (is_directory(it->status()))
    //         {
    //             std::string paramDir = it->path().string();
    //
    //             unsigned int dashpos = paramDir.find_last_of("/");
    //             std::string nameMesh = paramDir.substr(dashpos + 1);
    //
    //             if (M_meshASPairMap.find(nameMesh) == M_meshASPairMap.end())
    //             {
    //                 shp<TreeNode> defTreeNode = generateDefaultTreeNode(nameMesh);
    //                 shp<AssemblerType> defAssembler = AssemblerFactory<FEVECTOR COMMA FEMATRIX>(M_data, defTreeNode);
    //                 defAssembler->setup();
    //
    //                 M_meshASPairMap[nameMesh].first = defAssembler;
    //                 M_meshASPairMap[nameMesh].second.resize(defAssembler->getNumComponents());
    //
    //                 M_bases[nameMesh].reset(new RBBases(M_data, M_comm));
    //                 M_bases[nameMesh]->setNumberOfFields(defAssembler->getNumComponents());
    //                 unsigned int indexField = 0;
    //                 while (defAssembler->getFEspace(indexField))
    //                 {
    //                     M_bases[nameMesh]->setFESpace(defAssembler->getFEspace(indexField),
    //                                                   indexField);
    //                     indexField++;
    //                 }
    //             }
    //         }
    //     }
    //
    //     i++;
    // }
    // printlog(MAGENTA, "done\n", M_data.getVerbose());
}

void
BasisGenerator::
orthonormalize()
{
    // printlog(MAGENTA, "[BasisGenerator] orthonormalize basis ... \n", M_data.getVerbose());
    //
    // for (auto& meshas : M_meshASPairMap)
    // {
    //     auto curBasis = M_bases[meshas.first];
    //
    //     for (unsigned int i = 0; i < curBasis->getNumFields(); i++)
    //     {
    //         if (curBasis->hasSupremizers(i))
    //             curBasis->normalizeBasis(i, meshas.second.first->getNorm(i).data());
    //     }
    // }
}

void
BasisGenerator::
parseFiles()
{
    // printlog(MAGENTA, "[BasisGenerator] parsing files ... \n", M_data.getVerbose());
    //
    // // using namespace boost::filesystem;
    //
    // std::string snapshotsdir = M_data("rb/offline/snapshots/directory", "snapshots");
    //
    // if (!exists(snapshotsdir))
    //     throw new Exception("Snapshots directory has not been generated yet!");
    //
    // std::string paramdir = snapshotsdir + "/param";
    // unsigned int i = 0;
    // // we loop over the folders with the parameters
    // while (fs::exists(paramdir + std::to_string(i)))
    // {
    //     directory_iterator end_it;
    //     for (directory_iterator it(paramdir + std::to_string(i)); it != end_it; it++)
    //     {
    //         if (is_directory(it->status()))
    //             parseParameterSnapshots(it->path().string());
    //     }
    //
    //     i++;
    // }
    // printlog(MAGENTA, "done\n", M_data.getVerbose());
}

LifeV::LinearSolver
BasisGenerator::
setupLinearSolver(SparseMatrix matrix)
{

    // // solver part
    // LifeV::LinearSolver linearSolver(M_comm);
    // linearSolver.setOperator(matrix.data());
    //
    // Teuchos::RCP<Teuchos::ParameterList> aztecList =
    //                                    Teuchos::rcp(new Teuchos::ParameterList);
    //
    // std::string xmlfile = M_data("rb/offline/basis/xmlfile", "SolverParamList.xml");
    // aztecList = Teuchos::getParametersFromXmlFile(xmlfile);
    // linearSolver.setParameters(*aztecList);
    //
    // typedef LifeV::PreconditionerML         precML_type;
    // typedef shp<precML_type>    precMLPtr_type;
    // precML_type * precRawPtr;
    // precRawPtr = new precML_type;
    //
    // GetPot dummyDatafile;
    // precRawPtr->setDataFromGetPot(dummyDatafile, "precMLL");
    // shp<LifeV::Preconditioner> precPtr;
    // precPtr.reset(precRawPtr);
    //
    // linearSolver.setPreconditioner(precPtr);
    //
    // return linearSolver;
}

void
BasisGenerator::
addSupremizers()
{
    // std::string outdir = M_data("rb/offline/basis/directory", "basis");
    // bool addprimalsupremizers = M_data("rb/offline/basis/addprimalsupremizers", true);
    // if (addprimalsupremizers)
    // {
    //     printlog(MAGENTA, "[BasisGenerator] adding primal supremizers ... \n", M_data.getVerbose());
    //
    //     unsigned int field2augment = M_data("rb/offline/basis/primal_supremizers/field2augment", 0);
    //     unsigned int limitingfield = M_data("rb/offline/basis/primal_supremizers/limitingfield", 1);
    //
    //     for (auto& meshas : M_meshASPairMap)
    //     {
    //         printlog(GREEN, "mesh = " + meshas.first + "\n", M_data.getVerbose());
    //         // note: BCs. Here the matrix must have ones on the nodes corresponding
    //         // to dirichlet nodes. So, if norm = mass + stiffness (H1 norm),
    //         // than mass can have 1s on diagonal and stiffness 0s
    //         MatrixEp normMatrix = meshas.second.first->getNorm(field2augment);
    //
    //         // here the matrix must have 0 on the nodes corresponding to dirichlet
    //         // nodes
    //         MatrixEp constraintMatrix = meshas.second.first->getConstraintMatrix();
    //         constraintMatrix.dump(outdir + "/" + meshas.first + "/primalConstraint");
    //         auto map = *constraintMatrix.data()->rangeMapPtr();
    //         auto linearSolver = setupLinearSolver(normMatrix);
    //
    //         std::vector<shp<VECTOREPETRA>> basisFunctions = M_bases[meshas.first]->getFullBasis(limitingfield);
    //
    //         for (unsigned int i = 0; i < basisFunctions.size(); i++)
    //         {
    //             printlog(YELLOW, "adding supremizer " + std::to_string(i) + " ... \n",
    //                      M_data.getVerbose());
    //
    //             VectorEp basisFunction;
    //             basisFunction = basisFunctions[i];
    //
    //             VectorEp rhs = constraintMatrix * basisFunction;
    //
    //             linearSolver.setRightHandSide(rhs.data());
    //
    //             shp<VECTOREPETRA> solution(new VECTOREPETRA(map, LifeV::Unique));
    //             linearSolver.solve(solution);
    //
    //             M_bases[meshas.first]->addPrimalSupremizer(solution, field2augment, limitingfield);
    //         }
    //     }
    // }
    //
    // bool adddualsupremizers = M_data("rb/offline/basis/adddualsupremizers", true);
    // if (adddualsupremizers)
    // {
    //     printlog(MAGENTA, "[BasisGenerator] adding dual supremizers ... \n", M_data.getVerbose());
    //
    //     unsigned int field2augment = M_data("rb/offline/basis/dual_supremizers/field2augment", 0);
    //
    //     for (auto& meshas : M_meshASPairMap)
    //     {
    //         printlog(GREEN, "mesh = " + meshas.first + "\n", M_data.getVerbose());
    //         // note: BCs. same as above
    //         MatrixEp normMatrix = meshas.second.first->getNorm(field2augment);
    //
    //         InterfaceAssembler<VectorEp, MatrixEp> interfaceAssembler(M_data);
    //
    //         SHP(BuildingBlock) buildingBlock = meshas.second.first->getTreeNode()->M_block;
    //
    //         // create list of faces
    //         std::vector<GeometricFace> faces = buildingBlock->getOutlets();
    //         faces.push_back(buildingBlock->getInlet());
    //
    //         for (auto face : faces)
    //         {
    //             BlockMatrix<MatrixEp> constraintMatrixBlock;
    //             BlockMatrix<MatrixEp> constraintMatrixDummyBlock;
    //             interfaceAssembler.buildCouplingMatrices(meshas.second.first,
    //                                                      face,
    //                                                      constraintMatrixBlock,
    //                                                      constraintMatrixDummyBlock);
    //             // we assume that the first block is the one to be coupled
    //             // (as in interface assembler)
    //             MatrixEp constraintMatrix = constraintMatrixBlock.block(0,0);
    //             constraintMatrix.dump(outdir + "/" + meshas.first + "/dualConstraint" + std::to_string(face.M_flag));
    //             auto map = *constraintMatrix.data()->rangeMapPtr();
    //             auto linearSolver = setupLinearSolver(normMatrix);
    //             auto lagrangeMap = *constraintMatrix.data()->domainMapPtr();
    //             unsigned int numLagrange = lagrangeMap.mapSize();
    //
    //             for (unsigned int i = 0; i < numLagrange; i++)
    //             {
    //                 printlog(YELLOW, "adding supremizer " + std::to_string(i) + " ... \n",
    //                          M_data.getVerbose());
    //
    //                 shp<VECTOREPETRA> selector(new VECTOREPETRA(lagrangeMap, LifeV::Unique));
    //                 selector->zero();
    //                 selector->operator[](i) = 1.0;
    //
    //                 VectorEp basisFunction;
    //                 basisFunction = selector;
    //
    //                 VectorEp rhs = constraintMatrix * basisFunction;
    //
    //                 linearSolver.setRightHandSide(rhs.data());
    //
    //                 shp<VECTOREPETRA> solution(new VECTOREPETRA(map, LifeV::Unique));
    //                 linearSolver.solve(solution);
    //
    //                 M_bases[meshas.first]->addDualSupremizer(solution, field2augment);
    //             }
    //         }
    //     }
    // }
}

void
BasisGenerator::
performPOD()
{
    // using namespace rbLifeV;
    //
    // double podtol = M_data("rb/offline/basis/podtol", 1e-5);
    //
    // std::string outdir = M_data("rb/offline/basis/directory", "basis");
    // fs::create_directory(outdir);
    //
    // printlog(MAGENTA, "[BasisGenerator] performing POD ... \n", M_data.getVerbose());
    // for (auto pair : M_meshASPairMap)
    // {
    //     unsigned int count = 0;
    //     VectorFunctions newBasisFunctions(pair.second.second.size());
    //     fs::create_directory(outdir + "/" + pair.first);
    //
    //     for (auto sn : pair.second.second)
    //     {
    //         ProperOrthogonalDecomposition pod(M_comm,
    //                                           pair.second.first->getFEspace(count)->map(),
    //                                           true);
    //         pod.initPOD(sn.size(), sn.data(), pair.second.first->getNorm(count).data());
    //         pair.second.first->getNorm(count).dump(outdir + "/" + pair.first + "/norm" +
    //                            std::to_string(count));
    //         pod.setSvdFileName(outdir + "/" + pair.first + "/svd" +
    //                            std::to_string(count) + ".txt");
    //         pod.generatePODbasisTol(podtol);
    //
    //         unsigned int nbfs = pod.getRBdimension();
    //         std::vector<shp<VECTOREPETRA>> basisFunctions(nbfs);
    //
    //         pod.swapReducedBasis(basisFunctions, 0);
    //         M_bases[pair.first]->setPath(outdir + "/" + pair.first);
    //         M_bases[pair.first]->setBasisFunctions(basisFunctions, count);
    //         count++;
    //     }
    // }
    // printlog(MAGENTA, "done\n", M_data.getVerbose());
}

shp<TreeNode>
BasisGenerator::
generateDefaultTreeNode(const std::string& nameMesh)
{
    // if (nameMesh.find("tube") != std::string::npos)
    //     return generateDefaultTube(nameMesh);
    // else if (nameMesh.find("bifurcation_symmetric"))
    //     return generateDefaultSymmetricBifurcation(nameMesh);
    // else
    // {
    //     throw new Exception("BasisGenerator: this branch must still be implemented");
    // }
    //
    // return nullptr;
}

shp<TreeNode>
BasisGenerator::
generateDefaultTube(const std::string& nameMesh)
{
    // unsigned int diameter = std::atoi(nameMesh.substr(5,6).c_str());
    // unsigned int length = std::atoi(nameMesh.substr(7,8).c_str());
    // std::string refinement = nameMesh.substr(9);
    //
    // shp<Tube> defaultTube(new Tube(M_comm, refinement, false, diameter, length));
    // defaultTube->readMesh();
    //
    // shp<TreeNode> treeNode(new TreeNode(defaultTube, 1234));
    //
    // return treeNode;
}

shp<TreeNode>
BasisGenerator::
generateDefaultSymmetricBifurcation(const std::string& nameMesh)
{
    // unsigned int alpha = std::atoi(nameMesh.substr(13,15).c_str());
    // std::string refinement = nameMesh.substr(17);
    //
    // shp<BifurcationSymmetric> defaultBifurcation(new BifurcationSymmetric(M_comm, refinement, false, alpha));
    // defaultBifurcation->readMesh();
    //
    // shp<TreeNode> treeNode(new TreeNode(defaultBifurcation, 1234));
    //
    // return treeNode;
}

void
BasisGenerator::
parseParameterSnapshots(const std::string& paramDir)
{
    // // using namespace boost::filesystem;
    //
    // unsigned int dashpos = paramDir.find_last_of("/");
    // std::string nameMesh = paramDir.substr(dashpos + 1);
    //
    // directory_iterator end_it;
    // for (directory_iterator it(paramDir); it != end_it; it++)
    // {
    //     std::string curfile = it->path().string();
    //     if (curfile.find(".snap") != std::string::npos)
    //     {
    //         unsigned int dotpos = curfile.find_last_of(".");
    //         unsigned int componentIndex = std::atoi(curfile.substr(dotpos-1,dotpos).c_str());
    //
    //         addSnapshotsFromFile(curfile,
    //                              M_meshASPairMap[nameMesh].second[componentIndex],
    //                              M_meshASPairMap[nameMesh].first->getFEspace(componentIndex));
    //     }
    // }
}

void
BasisGenerator::
addSnapshotsFromFile(const std::string& snapshotsFile,
                     std::vector<shp<VECTOREPETRA>>& snapshots,
                     shp<FESPACE> fespace)
{

    // std::ifstream infile(snapshotsFile);
    // std::string line;
    // while(std::getline(infile,line))
    // {
    //     shp<VECTOREPETRA> newVector(new VECTOREPETRA(fespace->map()));
    //
    //     std::stringstream linestream(line);
    //     std::string value;
    //     unsigned int i = 0;
    //     while(getline(linestream,value,','))
    //     {
    //         newVector->operator[](i) = std::atof(value.c_str());
    //         i++;
    //     }
    //     if (i != newVector->epetraVector().GlobalLength())
    //         throw new Exception("Stored snapshot length does not match fespace dimension!");
    //
    //     snapshots.push_back(newVector);
    // }
    // infile.close();
}

void
BasisGenerator::
dumpBasis()
{
    // // using namespace boost::filesystem;
    //
    // std::string outdir = M_data("rb/offline/basis/directory", "basis");
    // bool binary = M_data("rb/offline/basis/dumpbinary", true);
    //
    // std::ios_base::openmode omode = std::ios_base::app;
    // if (binary)
    //     omode = omode | std::ios::binary;
    //
    // create_directory(outdir);
    //
    // for (auto meshbasis : M_bases)
    // {
    //     std::string meshdir = outdir + "/" + meshbasis.first;
    //     create_directory(meshdir);
    //
    //     M_bases[meshbasis.first]->dump();
    // }

}

}  // namespace RedMA
