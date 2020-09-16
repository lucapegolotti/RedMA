#include "SaddlePointPreconditioner.hpp"

namespace RedMA
{

SaddlePointPreconditioner::
SaddlePointPreconditioner(const DataContainer& data, const BM& matrix) :
  M_data(data),
  M_matrix(matrix)
{
    // Chrono chrono;
    // chrono.start();
    //
    // printlog(MAGENTA, "[SaddlePointPreconditioner] starting setup ...\n", M_data.getVerbose());
    // // split matrices
    // unsigned int nBlocks = matrix->nRows();
    //
    // // read options
    // setSolverOptions();
    //
    // // preconditioner used to approximate Am1 in iterations
    // M_innerPrecType = M_data("preconditioner/inner", "SIMPLE");
    //
    // // solution method to approximate BAm1BT in schur
    // M_approxSchurType = M_data("preconditioner/approxshur", "SIMPLE");
    //
    // if (nBlocks > 1)
    // {
    //     unsigned int nPrimal = 1;
    //     // if nBlocks == 2, then we are imposing weak boundary conditions and
    //     // we know that nprimal = 1
    //     if (nBlocks > 2)
    //     {
    //         while (nPrimal < nBlocks && matrix->block(0,nPrimal)->isZero())
    //             nPrimal++;
    //
    //         if (nPrimal == nBlocks)
    //             throw new Exception("The system has not a saddle point structure");
    //     }
    //     unsigned int nDual = nBlocks - nPrimal;
    //     M_nPrimalBlocks = nPrimal;
    //     M_nDualBlocks = nDual;
    //
    //     BM A = matrix->getSubmatrix(0, nPrimal-1, 0, nPrimal-1);
    //     BM BT = matrix->getSubmatrix(0, nPrimal-1, nPrimal, nBlocks-1);
    //     BM B = matrix->getSubmatrix(nPrimal, nBlocks-1, 0, nPrimal-1);
    //     BM C = matrix->getSubmatrix(nPrimal, nBlocks-1, nPrimal, nBlocks-1);
    //
    //     BlockMaps bmaps(*BT);
    //     M_primalMap = bmaps.getMonolithicRangeMapEpetra();
    //     M_dualMap = bmaps.getMonolithicDomainMapEpetra();
    //
    //     BlockMaps allmaps(*matrix);
    //     M_rangeMaps = allmaps.getRangeMapsEpetra();
    //     M_domainMaps = allmaps.getDomainMapsEpetra();
    //     M_monolithicMap = allmaps.getMonolithicRangeMapEpetra();
    //     M_setupTime = chrono.diff();
    //
    //     std::string msg = "Monolithic map size = ";
    //     msg += std::to_string(M_monolithicMap->mapSize());
    //     msg += "\n";
    //     printlog(GREEN, msg, M_data.getVerbose());
    //
    //     // this is to be optimized
    //     M_matrixCollapsed->softCopy(collapseBlocks(matrix, allmaps));
    //
    //     allocateInnerPreconditioners(A);
    //
    //     if (!std::strcmp(M_innerPrecType.c_str(), "exact") ||
    //         !std::strcmp(M_approxSchurType.c_str(), "exact"))
    //     {
    //         allocateInverseSolvers(A);
    //     }
    //
    //     printlog(MAGENTA, msg, M_data.getVerbose());
    //     computeSchurComplement(A, BT, B, C);
    //
    //     printlog(MAGENTA, msg, M_data.getVerbose());
    // }
    // else
    // {
    //     M_nPrimalBlocks = 1;
    //     M_nDualBlocks = 0;
    //     allocateInnerPreconditioners(matrix);
    // }
    //
    // M_setupTime = chrono.diff();
    //
    // std::string msg = "done, in ";
    // msg += std::to_string(M_setupTime);
    // msg += " seconds\n";
    // printlog(MAGENTA, msg, M_data.getVerbose());
}

void
SaddlePointPreconditioner::
computeSchurComplement(const BM& A, const BM& BT, const BM& B, const BM& C)
{
    // computeAm1BT(A, BT);
    //
    // M_S->softCopy(B->multiplyByMatrix(M_Am1BT));
    // M_S->multiplyByScalar(-1.0);
    // M_S->add(C);
    //
    // SHP(BlockMatrix) Scoll1;// = collapseBlocks(*M_S, BlockMaps(*M_S));
    // SHP(SparseMatrix) Scoll2;// = collapseBlocks(*Scoll1, BlockMaps(*Scoll1));
    //
    // // create invertible approximated matrix
    // M_approximatedSchurInverse.reset(new ApproxInv);
    //
    // SHP(Teuchos::ParameterList) globalSchurOptions;
    // globalSchurOptions.reset(new
    //     Teuchos::ParameterList(M_solversOptionsInner->sublist("GlobalSchurOperator")));
    //
    // M_approximatedSchurInverse->SetRowMatrix(static_cast<MATRIXEPETRA*>(Scoll2->data().get())->matrixPtr());
    // M_approximatedSchurInverse->SetParameterList(*globalSchurOptions);
    // M_approximatedSchurInverse->Compute();
}

void
SaddlePointPreconditioner::
allocateInverseSolvers(const BM& primalMatrix)
{
    // unsigned int nrows = primalMatrix->nRows();
    // for (unsigned int i = 0; i < nrows; i++)
    // {
    //     SHP(InvOp) invOper;
    //     SHP(NSOp) oper(new NSOp);
    //
    //     boost::numeric::ublas::matrix<SHP(MATRIXEPETRA::matrix_type)> matrixGrid(2,2);
    //     matrixGrid(0,0) = static_cast<MATRIXEPETRA*>(primalMatrix->block(i,i)->block(0,0)->data().get())->matrixPtr();
    //     matrixGrid(1,0) = static_cast<MATRIXEPETRA*>(primalMatrix->block(i,i)->block(1,0)->data().get())->matrixPtr();
    //     matrixGrid(0,1) = static_cast<MATRIXEPETRA*>(primalMatrix->block(i,i)->block(0,1)->data().get())->matrixPtr();
    //     if (!primalMatrix->block(i,i)->block(1,1)->isZero())
    //         matrixGrid(1,1) = static_cast<MATRIXEPETRA*>(primalMatrix->block(i,i)->block(1,1)->data().get())->matrixPtr();
    //
    //     oper->setUp(matrixGrid,
    //     static_cast<MATRIXEPETRA*>(primalMatrix->block(i,i)->block(0,0)->data().get())->rangeMap().commPtr());
    //
    //     SHP(Teuchos::ParameterList) options;
    //
    //     options.reset(new Teuchos::ParameterList(M_solversOptionsInner->sublist("InnerBlockOperator")));
    //
    //     double innertol = M_data("preconditioner/innertol", 1e-3);
    //
    //     std::string solverType(options->get<std::string>("Linear Solver Type"));
    //     if (!std::strcmp(solverType.c_str(),"AztecOO"))
    //         options->sublist(solverType).sublist("options")
    //                                     .get<double>("tol") = innertol;
    //     else if (!std::strcmp(solverType.c_str(),"Belos"))
    //         options->sublist(solverType).sublist("options")
    //                                     .get<double>("Convergence Tolerance") = innertol;
    //
    //     invOper.reset(LifeV::Operators::InvertibleOperatorFactory::instance().createObject(solverType));
    //
    //     invOper->setParameterList(options->sublist(solverType));
    //     invOper->setOperator(oper);
    //     invOper->setPreconditioner(M_innerPreconditioners[i]);
    //     M_invOperators.push_back(invOper);
    // }
}

void
SaddlePointPreconditioner::
allocateInnerPreconditioners(const BM& primalMatrix)
{
    // unsigned int nrows = primalMatrix->nRows();
    //
    // for (unsigned int i = 0; i < nrows; i++)
    // {
    //     SHP(NSPrec) newPrec;
    //     newPrec.reset(LifeV::Operators::NSPreconditionerFactory::
    //                   instance().createObject("SIMPLE"));
    //     newPrec->setOptions(*M_solversOptionsInner);
    //
    //     if (!primalMatrix->block(i,i)->block(1,1)->isZero())
    //     {
    //         newPrec->setUp(SHP(MATRIXEPETRA)(static_cast<MATRIXEPETRA*>(primalMatrix->block(i,i)->block(0,0)->data().get())),
    //                        SHP(MATRIXEPETRA)(static_cast<MATRIXEPETRA*>(primalMatrix->block(i,i)->block(1,0)->data().get())),
    //                        SHP(MATRIXEPETRA)(static_cast<MATRIXEPETRA*>(primalMatrix->block(i,i)->block(0,1)->data().get())),
    //                        SHP(MATRIXEPETRA)(static_cast<MATRIXEPETRA*>(primalMatrix->block(i,i)->block(1,1)->data().get())));
    //     }
    //     else
    //     {
    //         newPrec->setUp(SHP(MATRIXEPETRA)(static_cast<MATRIXEPETRA*>(primalMatrix->block(i,i)->block(0,0)->data().get())),
    //                        SHP(MATRIXEPETRA)(static_cast<MATRIXEPETRA*>(primalMatrix->block(i,i)->block(1,0)->data().get())),
    //                        SHP(MATRIXEPETRA)(static_cast<MATRIXEPETRA*>(primalMatrix->block(i,i)->block(0,1)->data().get())));
    //     }
    //
    //     std::vector<SHP(Epetra_Map)> localRangeMaps(2);
    //     std::vector<SHP(Epetra_Map)> localDomainMaps(2);
    //
    //     for (unsigned int j = 0; j < 2; j++)
    //         localRangeMaps[j].reset(new Epetra_Map(
    //         static_cast<MATRIXEPETRA*>(primalMatrix->block(i,i)->block(j,0)->data().get())->matrixPtr()->OperatorRangeMap()));
    //
    //     for (unsigned int j = 0; j < 2; j++)
    //         localDomainMaps[j].reset(new Epetra_Map(
    //         static_cast<MATRIXEPETRA*>(primalMatrix->block(i,i)->block(0,j)->data().get())->matrixPtr()->OperatorDomainMap()));
    //
    //     SHP(LifeV::BlockEpetra_Map) rmblock(new LifeV::BlockEpetra_Map(localRangeMaps));
    //     SHP(LifeV::BlockEpetra_Map) dmblock(new LifeV::BlockEpetra_Map(localDomainMaps));
    //     newPrec->setRangeMap(rmblock);
    //     newPrec->setDomainMap(dmblock);
    //
    //     newPrec->updateApproximatedMomentumOperator();
    //     newPrec->updateApproximatedSchurComplementOperator();
    //     M_innerPreconditioners.push_back(newPrec);
    // }
}

void
SaddlePointPreconditioner::
computeAm1BT(const BM& A, const BM& BT)
{
    // M_Am1BT->resize(A->nRows(), BT->nCols());
    //
    // for (unsigned int i = 0; i < A->nRows(); i++)
    // {
    //     for (unsigned int j = 0; j < BT->nCols(); j++)
    //     {
    //         BM Abm(dynamic_cast<BlockMatrix*>(A->block(i,i).get()));
    //         BM BTbm(dynamic_cast<BlockMatrix*>(BT->block(i,j).get()));
    //         M_Am1BT->block(i,j)->softCopy(computeSingleAm1BT(Abm, BTbm,i));
    //     }
    // }
    // M_Am1BT->close();
}

SHP(aMatrix)
SaddlePointPreconditioner::
computeSingleAm1BT(const BM& A, const BM& BT,
                   const unsigned int& index)
{
    // BlockMatrix retMat;
    // retMat.resize(A->nRows(), BT->nCols());
    //
    // if (!BT->isZero())
    // {
    //     Chrono chrono;
    //     chrono.start();
    //     printlog(GREEN, "[SaddlePointPreconditioner] single AM1BT ...", M_data.getVerbose());
    //
    //
    //     MAPEPETRA rangeMapU = static_cast<MATRIXEPETRA*>(BT->block(0,0)->data().get())->rangeMap();
    //     MAPEPETRA rangeMapP = static_cast<MATRIXEPETRA*>(A->block(1,0)->data().get())->rangeMap();
    //     MAPEPETRA domainMap = static_cast<MATRIXEPETRA*>(BT->block(0,0)->data().get())->domainMap();
    //
    //     VECTOREPETRA selector(domainMap);
    //     VECTOREPETRA colU(rangeMapU);
    //     VECTOREPETRA colP(rangeMapP);
    //
    //     unsigned int ncols = domainMap.mapSize();
    //
    //     std::vector<DistributedVector> ressU(ncols);
    //     std::vector<DistributedVector> ressP(ncols);
    //     for (unsigned int i = 0; i < ncols; i++)
    //     {
    //         selector.zero();
    //         colU.zero();
    //         colP.zero();
    //         if (domainMap.isOwned(i))
    //             selector[i] = 1.0;
    //
    //         colU = (*BT.block(0,0).data()) * selector;
    //
    //         ressU[i]->data().reset(new VECTOREPETRA(rangeMapU));
    //         ressP[i]->data().reset(new VECTOREPETRA(rangeMapP));
    //
    //         if (!std::strcmp(M_approxSchurType.c_str(),"exact"))
    //         {
    //             auto monolithicMap = rangeMapU;
    //             monolithicMap += rangeMapP;
    //
    //             VECTOREPETRA X(monolithicMap);
    //             VECTOREPETRA Y(monolithicMap);
    //             X.zero();
    //             Y.zero();
    //
    //             X.subset(colU, rangeMapU, 0, 0);
    //             CoutRedirecter ct;
    //             ct.redirect();
    //             M_invOperators[index]->ApplyInverse(X.epetraVector(), Y.epetraVector());
    //             ct.restore();
    //
    //             ressU[i].data()->subset(Y, rangeMapU, 0, 0);
    //             ressP[i].data()->subset(Y, rangeMapP, rangeMapU.mapSize(), 0);
    //         }
    //         else if (!std::strcmp(M_approxSchurType.c_str(),"SIMPLE"))
    //             M_innerPreconditioners[index]->ApplyInverse(colU, colP,
    //                                                         *ressU[i].data(),
    //                                                         *ressP[i].data());
    //
    //     }
    //     MatrixEp Am1BTu(ressU);
    //     MatrixEp Am1BTp(ressP);
    //     retMat.block(0,0).softCopy(Am1BTu);
    //     retMat.block(1,0).softCopy(Am1BTp);
    //
    //     std::string msg = "done, in ";
    //     msg += std::to_string(chrono.diff());
    //     msg += " seconds\n";
    //     printlog(GREEN, msg, M_data.getVerbose());
    // }
    // return retMat;
}

void
SaddlePointPreconditioner::
setSolverOptions()
{
    std::string optionsPrec = M_data("preconditioner/options",
                                     "datafiles/solversOptionsFast");
    optionsPrec += ".xml";
    M_solversOptionsInner = Teuchos::getParametersFromXmlFile(optionsPrec);
    std::shared_ptr<Teuchos::ParameterList> monolithicOptions;
    monolithicOptions.reset(
        new Teuchos::ParameterList(M_solversOptionsInner->sublist("InnerBlockOperator")));
    M_pListLinSolver = monolithicOptions;
}

void
SaddlePointPreconditioner::
solveEveryPrimalBlock(const VECTOREPETRA& X, VECTOREPETRA &Y) const
{
    // unsigned int offset = 0;
    // unsigned int count = 0;
    // for (unsigned int i = 0; i < M_nPrimalBlocks; i++)
    // {
    //     MAPEPETRA rangeVelocity = *M_rangeMaps[count];
    //     MAPEPETRA rangePressure = *M_rangeMaps[count+1];
    //
    //     if (!std::strcmp(M_innerPrecType.c_str(), "exact"))
    //     {
    //         auto monolithicMap = rangeVelocity;
    //         monolithicMap += rangePressure;
    //
    //         VECTOREPETRA subX(monolithicMap);
    //         VECTOREPETRA subY(monolithicMap);
    //
    //         subX.subset(X, monolithicMap, offset, 0);
    //
    //         RedMA::CoutRedirecter ct;
    //         ct.redirect();
    //
    //         M_invOperators[i]->ApplyInverse(subX.epetraVector(),
    //                                         subY.epetraVector());
    //
    //         ct.restore();
    //
    //         Y.subset(subY, monolithicMap, 0, offset);
    //     }
    //     else if (!std::strcmp(M_innerPrecType.c_str(), "SIMPLE"))
    //     {
    //         VECTOREPETRA subVelocity(rangeVelocity);
    //         VECTOREPETRA subPressure(rangePressure);
    //         subVelocity.zero();
    //         subPressure.zero();
    //
    //         subVelocity.subset(X, rangeVelocity, offset, 0);
    //         subPressure.subset(X, rangePressure,
    //                            offset + rangeVelocity.mapSize(), 0);
    //
    //         VECTOREPETRA resVelocity(rangeVelocity);
    //         VECTOREPETRA resPressure(rangePressure);
    //
    //         M_innerPreconditioners[i]->ApplyInverse(subVelocity, subPressure,
    //                                                 resVelocity, resPressure);
    //
    //         Y.subset(resVelocity, rangeVelocity, 0, offset);
    //         Y.subset(resPressure, rangePressure, 0,
    //                  offset + rangeVelocity.mapSize());
    //     }
    //     offset += rangeVelocity.mapSize() + rangePressure.mapSize();
    //     count = count + 2;
    // }
}

void
SaddlePointPreconditioner::
applyEveryB(const VECTOREPETRA& X, VECTOREPETRA &Y) const
{
    // unsigned int offset = 0;
    //
    // std::vector<VECTOREPETRA> Xs;
    // for (unsigned int i = 0; i < M_nPrimalBlocks * 2; i++)
    // {
    //     VECTOREPETRA newVec(*M_rangeMaps[i]);
    //     newVec.subset(X, *M_rangeMaps[i], offset, 0);
    //     Xs.push_back(newVec);
    //     offset += M_rangeMaps[i]->mapSize();
    // }
    //
    // offset = 0;
    // for (unsigned int i = M_nPrimalBlocks * 2; i < M_nPrimalBlocks * 2 + M_nDualBlocks; i++)
    // {
    //     MAPEPETRA curRange = *M_rangeMaps[i];
    //     VECTOREPETRA subRes(curRange, LifeV::Unique);
    //     subRes.zero();
    //     for (unsigned int j = 0; j < M_nPrimalBlocks * 2; j++)
    //     {
    //         if (!M_matrixCollapsed.block(i,j).isNull())
    //             subRes += (*M_matrixCollapsed.block(i,j).data()) * Xs[j];
    //     }
    //     Y.subset(subRes, curRange, 0, offset);
    //     offset += curRange.mapSize();
    // }
}

void
SaddlePointPreconditioner::
applyEveryAm1BT(const VECTOREPETRA& X, VECTOREPETRA &Y) const
{
    // unsigned int offset = 0;
    //
    // std::vector<VECTOREPETRA> Xs;
    // for (unsigned int i = M_nPrimalBlocks * 2; i < M_nPrimalBlocks * 2 + M_nDualBlocks; i++)
    // {
    //     VECTOREPETRA newVec(*M_domainMaps[i]);
    //     newVec.subset(X, *M_domainMaps[i], offset, 0);
    //     Xs.push_back(newVec);
    //     offset += M_domainMaps[i]->mapSize();
    // }
    //
    // offset = 0;
    // for (unsigned int i = 0; i < 2 * M_nPrimalBlocks; i++)
    // {
    //     MAPEPETRA curRange = *M_domainMaps[i];
    //     VECTOREPETRA subRes(curRange, LifeV::Unique);
    //     subRes.zero();
    //     // for (unsigned int j = 2 * M_nPrimalBlocks; j < M_nPrimalBlocks * 2 + M_nDualBlocks; j++)
    //     // {
    //     //     if (!M_matrixCollapsed.block(i,j).isNull())
    //     //     {
    //     //         subRes += (*M_matrixCollapsed.block(i,j).data()) * Xs[j-2*M_nPrimalBlocks];
    //     //     }
    //     // }
    //     for (unsigned int j = 0; j < M_nDualBlocks; j++)
    //     {
    //         if (!M_Am1BT.block(i/2,j).block(i%2,0).isNull())
    //         {
    //             subRes += (*M_Am1BT.block(i/2,j).block(i%2,0).data()) * Xs[j];
    //         }
    //     }
    //     Y.subset(subRes, curRange, 0, offset);
    //     offset += curRange.mapSize();
    // }
}


int
SaddlePointPreconditioner::
ApplyInverse(const super::vector_Type& X, super::vector_Type& Y) const
{
    // if (M_nDualBlocks > 0)
    // {
    //     const VECTOREPETRA X_vectorEpetra(X, M_monolithicMap, LifeV::Unique);
    //     VECTOREPETRA Y_vectorEpetra(M_monolithicMap, LifeV::Unique);
    //
    //     VECTOREPETRA X_primal(M_primalMap, LifeV::Unique);
    //     VECTOREPETRA X_dual(M_dualMap, LifeV::Unique);
    //     X_primal.subset(X_vectorEpetra, *M_primalMap, 0, 0);
    //     X_dual.subset(X_vectorEpetra, *M_dualMap, M_primalMap->mapSize(), 0);
    //
    //     // here we store the result
    //     VECTOREPETRA Y_primal(M_primalMap, LifeV::Unique);
    //     VECTOREPETRA Y_dual(M_dualMap, LifeV::Unique);
    //
    //     VECTOREPETRA Z(M_primalMap, LifeV::Unique);
    //     solveEveryPrimalBlock(X_primal, Z);
    //
    //
    //     VECTOREPETRA Bz(M_dualMap, LifeV::Unique);
    //     applyEveryB(Z, Bz);
    //
    //     M_approximatedSchurInverse->ApplyInverse((X_dual-Bz).epetraVector(),
    //                                              Y_dual.epetraVector());
    //
    //
    //     VECTOREPETRA Am1BTy(M_primalMap, LifeV::Unique);
    //     applyEveryAm1BT(Y_dual, Am1BTy);
    //
    //     // this can be optimized because we already have Am1BTy
    //     // solveEveryPrimalBlock(BTy, Am1BTy);
    //
    //     Y_primal = Z - Am1BTy;
    //     Y_vectorEpetra.subset(Y_primal, *M_primalMap, 0, 0);
    //     Y_vectorEpetra.subset(Y_dual, *M_dualMap, 0, M_primalMap->mapSize());
    //     Y = dynamic_cast<Epetra_MultiVector&>(Y_vectorEpetra.epetraVector());
    // }
    // else
    // {
    //     M_innerPreconditioners[0]->ApplyInverse(X, Y);
    // }

    return 0;
}

}
