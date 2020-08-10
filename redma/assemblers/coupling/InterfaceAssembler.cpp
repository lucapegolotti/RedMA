#include "InterfaceAssembler.hpp"

namespace RedMA
{

SHP(LifeV::QuadratureRule)
InterfaceAssembler::
generateQuadratureRule(std::string tag) const
{
  //   using namespace LifeV;
  //
  // //   0.33333333333333331       0.33333333333333331
  // // 2.06349616025259287E-002  0.48968251919873701
  // // 0.48968251919873701       2.06349616025259287E-002
  // // 0.48968251919873701       0.48968251919873701
  // // 0.12582081701412900       0.43708959149293553
  // // 0.43708959149293553       0.12582081701412900
  // // 0.43708959149293553       0.43708959149293553
  // // 0.62359292876193562       0.18820353561903219
  // // 0.18820353561903219       0.62359292876193562
  // // 0.18820353561903219       0.18820353561903219
  // // 0.91054097321109406       4.47295133944529688E-002
  // // 4.47295133944529688E-002  0.91054097321109406
  // // 4.47295133944529688E-002  4.47295133944529688E-002
  // // 0.74119859878449801       3.68384120547362581E-002
  // // 0.74119859878449801       0.22196298916076573
  // // 3.68384120547362581E-002  0.74119859878449801
  // // 3.68384120547362581E-002  0.22196298916076573
  // // 0.22196298916076573       0.74119859878449801
  // // 0.22196298916076573       3.68384120547362581E-002
  //
  // // 9.71357962827961025E-002
  // // 3.13347002271398278E-002
  // // 3.13347002271398278E-002
  // // 3.13347002271398278E-002
  // // 7.78275410047754301E-002
  // // 7.78275410047754301E-002
  // // 7.78275410047754301E-002
  // // 7.96477389272090969E-002
  // // 7.96477389272090969E-002
  // // 7.96477389272090969E-002
  // // 2.55776756586981006E-002
  // // 2.55776756586981006E-002
  // // 2.55776756586981006E-002
  // // 4.32835393772893970E-002
  // // 4.32835393772893970E-002
  // // 4.32835393772893970E-002
  // // 4.32835393772893970E-002
  // // 4.32835393772893970E-002
  // // 4.32835393772893970E-002
  //
  //   SHP(QuadratureRule) customRule;
  //   if (!std::strcmp(tag.c_str(),"STRANG10"))
  //   {
  //       // rule taken from https://people.sc.fsu.edu/~jburkardt/datasets/quadrature_rules_tri/quadrature_rules_tri.html
  //       customRule.reset(new QuadratureRule("STRANG10",TRIANGLE, 3, 7, 0));
  //
  //       QuadraturePoint p1(0.333333333333333,0.333333333333333,0,
  //                          -0.149570044467670/2);
  //       customRule->addPoint(p1);
  //
  //       QuadraturePoint p2(0.479308067841923,0.260345966079038,0,
  //                          0.175615257433204/2);
  //       customRule->addPoint(p2);
  //
  //       QuadraturePoint p3(0.260345966079038,0.479308067841923,0,
  //                          0.175615257433204/2);
  //       customRule->addPoint(p3);
  //
  //       QuadraturePoint p4(0.260345966079038,0.260345966079038,0,
  //                          0.175615257433204/2);
  //       customRule->addPoint(p4);
  //
  //       QuadraturePoint p5(0.869739794195568,0.065130102902216,0,
  //                          0.053347235608839/2);
  //       customRule->addPoint(p5);
  //
  //       QuadraturePoint p6(0.065130102902216,0.869739794195568,0,
  //                          0.053347235608839/2);
  //       customRule->addPoint(p6);
  //
  //       QuadraturePoint p7(0.065130102902216,0.065130102902216,0,
  //                          0.053347235608839/2);
  //       customRule->addPoint(p7);
  //
  //       QuadraturePoint p8(0.638444188569809,0.312865496004875,0,
  //                          0.077113760890257/2);
  //       customRule->addPoint(p8);
  //
  //       QuadraturePoint p9(0.638444188569809,0.048690315425316,0,
  //                          0.077113760890257/2);
  //       customRule->addPoint(p9);
  //
  //       QuadraturePoint p10(0.312865496004875,0.638444188569809,0,
  //                          0.077113760890257/2);
  //       customRule->addPoint(p10);
  //
  //       QuadraturePoint p11(0.312865496004875,0.048690315425316,0,
  //                          0.077113760890257/2);
  //       customRule->addPoint(p11);
  //
  //       QuadraturePoint p12(0.048690315425316,0.638444188569809,0,
  //                          0.077113760890257/2);
  //       customRule->addPoint(p12);
  //
  //       QuadraturePoint p13(0.048690315425316,0.312865496004875,0,
  //                          0.077113760890257/2);
  //       customRule->addPoint(p13);
  //   }
  //   else if (!std::strcmp(tag.c_str(),"TOMS612_19"))
  //   {
  //       // rule taken from https://people.sc.fsu.edu/~jburkardt/datasets/quadrature_rules_tri/quadrature_rules_tri.html
  //       customRule.reset(new QuadratureRule("TOMS612_19",TRIANGLE, 3, 9, 0));
  //
  //       QuadraturePoint p1(0.33333333333333331,0.33333333333333331,0,
  //                          9.71357962827961025e-2/2);
  //       customRule->addPoint(p1);
  //
  //       QuadraturePoint p2(2.06349616025259287e-2,0.48968251919873701,0,
  //                          3.13347002271398278e-2/2);
  //       customRule->addPoint(p2);
  //
  //       QuadraturePoint p3(0.48968251919873701,2.06349616025259287e-2,0,
  //                          3.13347002271398278e-2/2);
  //       customRule->addPoint(p3);
  //
  //       QuadraturePoint p4(0.48968251919873701,0.48968251919873701,0,
  //                          3.13347002271398278e-2/2);
  //       customRule->addPoint(p4);
  //
  //       QuadraturePoint p5(0.12582081701412900,0.43708959149293553,0,
  //                          7.78275410047754301e-2/2);
  //       customRule->addPoint(p5);
  //
  //       QuadraturePoint p6(0.43708959149293553,0.12582081701412900,0,
  //                          7.78275410047754301e-2/2);
  //       customRule->addPoint(p6);
  //
  //       QuadraturePoint p7(0.43708959149293553,0.43708959149293553,0,
  //                          7.78275410047754301e-2/2);
  //       customRule->addPoint(p7);
  //
  //       QuadraturePoint p8(0.62359292876193562,0.18820353561903219,0,
  //                          7.96477389272090969e-2/2);
  //       customRule->addPoint(p8);
  //
  //       QuadraturePoint p9(0.18820353561903219,0.62359292876193562,0,
  //                          7.96477389272090969e-2/2);
  //       customRule->addPoint(p9);
  //
  //       QuadraturePoint p10(0.18820353561903219,0.18820353561903219,0,
  //                          7.96477389272090969e-2/2);
  //       customRule->addPoint(p10);
  //
  //       QuadraturePoint p11(0.91054097321109406,4.47295133944529688e-2,0,
  //                          2.55776756586981006e-2/2);
  //       customRule->addPoint(p11);
  //
  //       QuadraturePoint p12(4.47295133944529688e-2,0.91054097321109406,0,
  //                          2.55776756586981006e-2/2);
  //       customRule->addPoint(p12);
  //
  //       QuadraturePoint p13(4.47295133944529688e-2,4.47295133944529688e-2,0,
  //                          2.55776756586981006e-2/2);
  //       customRule->addPoint(p13);
  //
  //       //
  //
  //       QuadraturePoint p14(0.74119859878449801,3.68384120547362581e-2,0,
  //                          4.32835393772893970e-2/2);
  //       customRule->addPoint(p14);
  //
  //       QuadraturePoint p15(0.74119859878449801,0.22196298916076573,0,
  //                          4.32835393772893970e-2/2);
  //       customRule->addPoint(p15);
  //
  //       QuadraturePoint p16(3.6838412054736258e-2,0.74119859878449801,0,
  //                          4.32835393772893970e-2/2);
  //       customRule->addPoint(p1);
  //
  //       QuadraturePoint p17(3.68384120547362581e-2,0.22196298916076573,0,
  //                          4.32835393772893970e-2/2);
  //       customRule->addPoint(p17);
  //
  //       QuadraturePoint p18(0.22196298916076573,0.74119859878449801,0,
  //                          4.32835393772893970e-2/2);
  //       customRule->addPoint(p18);
  //
  //       QuadraturePoint p19(0.22196298916076573,3.68384120547362581e-2,0,
  //                          4.32835393772893970e-2/2);
  //       customRule->addPoint(p19);
  //   }
  //   else
  //   {
  //       throw new Exception("Quadrature rule " + tag + " not implemented!");
  //   }
  //   return customRule;
}

Interface::
Interface()
{
}

Interface::
Interface(SHP(AssemblerType) assemblerFather, const int& indexFather,
          SHP(AssemblerType) assemblerChild, const int& indexChild,
          const unsigned int& interfaceID) :
  M_assemblerFather(assemblerFather),
  M_indexFather(indexFather),
  M_assemblerChild(assemblerChild),
  M_indexChild(indexChild),
  M_ID(interfaceID)
{

}

InterfaceAssembler::
InterfaceAssembler(const DataContainer& data) :
  M_data(data),
  M_isInlet(false)
{
}

InterfaceAssembler::
InterfaceAssembler(const DataContainer& data,
                   const Interface& interface) :
  M_data(data),
  M_interface(interface),
  M_isInlet(false)
{
    // if (M_interface.M_indexFather == -1 || M_interface.M_indexChild == -1)
    //     M_isInlet = true;
    // setup();
}

void
InterfaceAssembler::
setup()
{
    // Chrono chrono;
    // chrono.start();
    //
    // printlog(YELLOW, "[InterfaceAssembler] initialize interface"
    //                  " assembler ...", M_data.getVerbose());
    //
    // M_stabilizationCoupling = M_data("coupling/stab_coefficient", 0.0);
    //
    // buildCouplingMatrices();
    //
    // std::string msg = "done, in ";
    // msg += std::to_string(chrono.diff());
    // msg += " seconds\n";
    // printlog(YELLOW, msg, M_data.getVerbose());
}

std::vector<SHP(aVector)>
InterfaceAssembler::
buildCouplingVectors(SHP(BasisFunctionFunctor) bfs,
                     const GeometricFace& face,
                     SHP(aAssembler) assembler) const
{
    // using namespace LifeV;
    // using namespace ExpressionAssembly;
    //
    // QuadratureBoundary boundaryQuadRule(buildTetraBDQR
    //                                    (*generateQuadratureRule("STRANG10")));
    //
    // // QuadratureBoundary boundaryQuadRule(buildTetraBDQR
    // //                                     (*generateQuadratureRule("TOMS612_19")));
    // // QuadratureBoundary boundaryQuadRule (buildTetraBDQR(quadRuleTria7pt));
    //
    // unsigned int nBasisFunctions = bfs->getNumBasisFunctions();
    //
    // std::vector<VectorEp> couplingVectors(3 * nBasisFunctions);
    //
    // SHP(ETFESPACE3) etfespace = assembler->getETFESpaceCoupling();
    // MAPEPETRA map = etfespace->map();
    // SHP(MESH) mesh = etfespace->mesh();
    // unsigned int faceFlag = face.M_flag;
    //
    // unsigned int count = 0;
    // // this must be changed for scalar problems (e.g. laplacian)
    // for (unsigned int dim = 0; dim < 3; dim++)
    // {
    //     LifeV::VectorSmall<3> versor;
    //     versor[0] = 0.;
    //     versor[1] = 0.;
    //     versor[2] = 0.;
    //
    //     versor[dim] = 1.;
    //     for (unsigned int i = 0; i < nBasisFunctions; i++)
    //     {
    //         SHP(VECTOREPETRA) currentMode(new VECTOREPETRA(map, LifeV::Repeated));
    //
    //         bfs->setIndex(i);
    //         integrate(boundary(mesh, faceFlag),
    //                   boundaryQuadRule,
    //                   etfespace,
    //                   eval(bfs, X) * dot(versor, phi_i)
    //               ) >> currentMode;
    //         couplingVectors[count].data() = currentMode;
    //         count++;
    //     }
    // }
    //
    // return couplingVectors;
}

void
InterfaceAssembler::
buildCouplingMatrices(SHP(AssemblerType) assembler, const GeometricFace& face,
                      BlockMatrix& matrixT, BlockMatrix& matrix)
{
    // SHP(BasisFunctionFunctor) bfs;
    //
    // bfs =  BasisFunctionFactory(M_data.getDatafile(), face, M_isInlet);
    //
    // std::vector<VectorEp> couplingVectors;
    // couplingVectors = buildCouplingVectors(bfs, face, assembler);
    //
    // matrixT.resize(assembler->getNumComponents(), 1);
    // matrix.resize(1, assembler->getNumComponents());
    // // we assume that the first field is the one to be coupled
    // matrixT.block(0,0).softCopy(MatrixEp(couplingVectors));
    // matrix.block(0,0).softCopy(matrixT.block(0,0).transpose());
    //
    // assembler->getBCManager()->apply0DirichletMatrix(matrixT,
    //                                                  assembler->getFESpaceBCs(),
    //                                                  assembler->getComponentBCs(),
    //                                                  0.0);
}

void
InterfaceAssembler::
addContributionRhs(const double& time, BlockVector& rhs,
                   const BlockVector& sol, const unsigned int& nPrimalBlocks)
{
    // unsigned int fatherID = M_interface.M_indexFather;
    // unsigned int childID = M_interface.M_indexChild;
    // unsigned int interfaceID = M_interface.M_ID;
    // SHP(aAssembler<InVectorType COMMA InMatrixType>) assemblerFather;
    // SHP(aAssembler<InVectorType COMMA InMatrixType>) assemblerChild;
    // assemblerFather = M_interface.M_assemblerFather;
    // assemblerChild = M_interface.M_assemblerChild;
    //
    // // we have (-1) because we are solving H un+1 = F(.) and coupling is in F
    // rhs.block(fatherID) -= M_fatherBT * sol.block(nPrimalBlocks + interfaceID);
    // rhs.block(childID)  -= M_childBT * sol.block(nPrimalBlocks + interfaceID);
    //
    // // no need to apply bcs as matrices have already bcs in them
    // // assemblerFather->getBCManager()->apply0DirichletBCs(rhs.block(fatherID),
    // //                                                     assemblerFather->getFESpaceBCs(),
    // //                                                     assemblerFather->getComponentBCs());
    // //
    // // assemblerChild->getBCManager()->apply0DirichletBCs(rhs.block(childID),
    // //                                                    assemblerChild->getFESpaceBCs(),
    // //                                                    assemblerChild->getComponentBCs());
    //
    // rhs.block(nPrimalBlocks + interfaceID) -= M_fatherB * sol.block(fatherID);
    // rhs.block(nPrimalBlocks + interfaceID) -= M_childB * sol.block(childID);
    //
    // if (M_stabilizationCoupling > THRESHOLDSTAB)
    // {
    //     rhs.block(nPrimalBlocks + interfaceID) -= (M_stabFather * sol.block(fatherID)) * (0.5 * M_stabilizationCoupling);
    //     rhs.block(nPrimalBlocks + interfaceID) -= (M_stabChild * sol.block(childID)) * (0.5 * M_stabilizationCoupling);
    //
    //     rhs.block(nPrimalBlocks + interfaceID) -=
    //     sol.block(nPrimalBlocks + interfaceID) * M_stabilizationCoupling;
    // }
}

double
InterfaceAssembler::
checkStabilizationTerm(const BlockVector& sol, const unsigned int& nPrimalBlocks)
{
    // if (M_stabilizationCoupling > THRESHOLDSTAB &&
    //     M_interface.M_assemblerFather && M_interface.M_assemblerChild)
    // {
    //     unsigned int fatherID = M_interface.M_indexFather;
    //     unsigned int childID = M_interface.M_indexChild;
    //     unsigned int interfaceID = M_interface.M_ID;
    //
    //     BlockVector<BlockVector<InVectorType>> res;
    //     res.resize(1);
    //
    //     res.block(0) -= (M_stabFather * sol.block(fatherID)) * 0.5;
    //     res.block(0) -= (M_stabChild * sol.block(childID)) * 0.5;
    //
    //     // std::cout << "--------" << std::endl << std::flush;
    //     // res.block(0).block(0).data()->showMe();
    //     // std::cout << "++++++++" << std::endl << std::flush;
    //     // sol.block(nPrimalBlocks + interfaceID).block(0).data()->showMe();
    //     //
    //     // std::cout << "stab term stress " << res.norm2() << std::endl << std::flush;
    //     // std::cout << "stab term lagrange " << sol.block(nPrimalBlocks + interfaceID).norm2() << std::endl << std::flush;
    //
    //     res.block(0) -= sol.block(nPrimalBlocks + interfaceID);
    //
    //     std::string msg = "[InterfaceAssembler] interface ID = ";
    //     msg += std::to_string(interfaceID);
    //     msg += ", stab term norm = ";
    //     msg += std::to_string(res.norm2());
    //     msg += "\n";
    //     printlog(MAGENTA, msg, M_data.getVerbose());
    // }
}

void
InterfaceAssembler::
addContributionJacobianRhs(const double& time, BlockMatrix& jac,
                           const BlockVector& sol, const unsigned int& nPrimalBlocks)
{
    // unsigned int fatherID = M_interface.M_indexFather;
    // unsigned int childID = M_interface.M_indexChild;
    // unsigned int interfaceID = M_interface.M_ID;
    //
    // // hard copy, otherwise we flip the sign of the matrices every time this
    // // function is called
    // jac.block(fatherID, nPrimalBlocks + interfaceID).hardCopy(M_fatherBT);
    // jac.block(childID,  nPrimalBlocks + interfaceID).hardCopy(M_childBT);
    // jac.block(nPrimalBlocks + interfaceID, fatherID).hardCopy(M_fatherB);
    // jac.block(nPrimalBlocks + interfaceID,  childID).hardCopy(M_childB);
    //
    // jac.block(fatherID, nPrimalBlocks + interfaceID) *= (-1);
    // jac.block(childID,  nPrimalBlocks + interfaceID) *= (-1);
    // jac.block(nPrimalBlocks + interfaceID, fatherID) *= (-1);
    // jac.block(nPrimalBlocks + interfaceID,  childID) *= (-1);
    //
    // if (M_stabilizationCoupling > THRESHOLDSTAB)
    // {
    //     jac.block(nPrimalBlocks + interfaceID, fatherID) += (M_stabFather * (-0.5 * M_stabilizationCoupling));
    //     jac.block(nPrimalBlocks + interfaceID,  childID) += (M_stabChild * (-0.5 * M_stabilizationCoupling));
    //
    //     jac.block(nPrimalBlocks + interfaceID, nPrimalBlocks + interfaceID).hardCopy(M_identity * (-1.0 * M_stabilizationCoupling));
    // }
}

}
