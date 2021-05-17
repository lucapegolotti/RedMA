#include "LinearSystemSolver.hpp"

namespace RedMA
{

LinearSystemSolver::
LinearSystemSolver(const DataContainer& data) :
  M_data(data),
  M_numSolves(0)
{
}

void
LinearSystemSolver::
solve(const BM& matrix, const BV& rhs, BV& sol)
{
    if (matrix->type() != BLOCK || rhs->type() != BLOCK)
        throw new Exception("LinearSystemSolver currently supports only block structures!");

    if (convert<BlockMatrix>(matrix)->block(0,0)->type() == DOUBLE)
    {
        double matValue = convert<Double>(convert<BlockMatrix>(matrix)->block(0,0))->getValue();
        double rhsValue = convert<Double>(convert<BlockVector>(rhs)->block(0))->getValue();
        sol->deepCopy(rhs);
        spcast<Double>(convert<BlockVector>(sol)->block(0))->setValue(rhsValue/matValue);
        return;
    }

    // CURRENTLY NOT WORKING!!
    /*if (convert<BlockMatrix>(matrix)->globalTypeIs(DENSE) &&
        convert<BlockVector>(rhs)->globalTypeIs(DENSE))
    {
        solveDense(matrix, rhs, sol);
        return;
    }*/

    BM matrixSparse = matrix;
    if (!convert<BlockMatrix>(matrixSparse)->globalTypeIs(SPARSE))
        matrixSparse = convert<BlockMatrix>(matrixSparse)->convertInnerTo(SPARSE,M_comm);

    BV rhsSparse = rhs;
    if(!convert<BlockVector>(rhsSparse)->globalTypeIs(DISTRIBUTED))
        rhsSparse = convert<BlockVector>(rhsSparse)->convertInnerTo(DISTRIBUTED, M_comm);

    //  blockMatrixToSparseMatrix(convert<BlockMatrix>(matrixSparse))->dump("matrixFEM");
    // blockVectorToDenseVector(convert<BlockVector>(rhsSparse))->dump("rhsFEM");

    if (!M_maps)
        M_maps.reset(new BlockMaps(convert<BlockMatrix>(matrixSparse)));
    else
        M_maps->updateCollapsedMatrix(convert<BlockMatrix>(matrixSparse));

    M_oper.reset(new LinearOperator(matrixSparse,M_maps));

    M_invOper.reset(new InverseOperator(M_data));
    M_invOper->setOperator(M_oper);
    M_invOper->setBlockMaps(M_maps);
    buildPreconditioner(matrixSparse);
    M_invOper->setPreconditioner(M_prec);

    Chrono chrono;
    chrono.start();
    printlog(MAGENTA, "[LinearSystemSolver] solve ...", M_data.getVerbose());

    M_statistics.M_numIterations = M_invOper->invert(rhsSparse, sol);

    M_statistics.M_solveTime = chrono.diff();
    std::string msg = "done, in ";
    msg += std::to_string(M_statistics.M_solveTime);
    msg += " seconds\n";
    printlog(GREEN, msg, M_data.getVerbose());

    M_numSolves++;
}

void
LinearSystemSolver::
buildPreconditioner(const BM& matrix)
{
    std::string precType = M_data("preconditioner/outer", "saddlepoint");

    if (!std::strcmp(precType.c_str(), "saddlepoint"))
    {
        unsigned int recomputeevery = M_data("preconditioner/recomputeevery", 1);
        if (M_prec == nullptr || (M_numSolves % recomputeevery) == 0)
            M_prec.reset(new SaddlePointPreconditioner(M_data,
                                                       spcast<BlockMatrix>(matrix),
                                                       spcast<BlockMatrix>(M_Mp)));
        else
            spcast<SaddlePointPreconditioner>(M_prec)->setup(spcast<BlockMatrix>(matrix),
                                                                 spcast<BlockMatrix>(M_Mp),
                                                                 false);
    }
    else
    {
        throw new Exception("Unrecognized type of preconditioner " + precType + " !");
    }
    M_statistics.M_precSetupTime = M_prec->getSetupTime();
}

void
LinearSystemSolver::
computeSchurComplementDense(const BM& matrix)
{
    shp<BlockMatrix> blockMatrix = convert<BlockMatrix>(matrix);

    // count primal blocks
    unsigned int nPrimal = 1;
    unsigned int nBlocks = matrix->nRows();

    while (blockMatrix->block(0, nPrimal)->isZero())
     nPrimal++;

    shp<BlockMatrix> A = convert<BlockMatrix>(blockMatrix->getSubmatrix(0, nPrimal-1, 0, nPrimal-1));
    shp<BlockMatrix> BT = convert<BlockMatrix>(blockMatrix->getSubmatrix(0, nPrimal-1, nPrimal, nBlocks-1));
    shp<BlockMatrix> B = convert<BlockMatrix>(blockMatrix->getSubmatrix(nPrimal, nBlocks-1, 0, nPrimal-1));
    shp<BlockMatrix> C = convert<BlockMatrix>(blockMatrix->getSubmatrix(nPrimal, nBlocks-1, nPrimal, nBlocks-1));

    M_collapsedAs.resize(nPrimal);
    M_solversAs.resize(nPrimal);

    for (unsigned int i = 0; i < nPrimal; i++)
    {
     M_solversAs[i].reset(new Epetra_SerialDenseSolver());
     M_collapsedAs[i] = blockMatrixToDenseMatrix(convert<BlockMatrix>(A->block(i,i)));
     M_solversAs[i]->SetMatrix(*spcast<DENSEMATRIX>(M_collapsedAs[i]->data()));
     M_solversAs[i]->Factor();
    }

    // compute Am1BT
    shp<BlockMatrix> Am1BT(new BlockMatrix());
    Am1BT->resize(BT->nRows(), BT->nCols());  // 1x1
    for (unsigned int i = 0; i < BT->nRows(); i++)
    {
        unsigned int v_dim = A->block(i,i)->block(0,0)->nRows();
        unsigned int p_dim = A->block(i,i)->block(1,0)->nRows();

        auto collapsedA = blockMatrixToDenseMatrix(convert<BlockMatrix>(A->block(i, i)));
        for (unsigned int j = 0; j < BT->nCols(); j++)
        {
            unsigned int nsrows = BT->block(i, j)->nRows();
            unsigned int nscols = BT->block(i, j)->nCols();

            shp<BlockMatrix> singleAm1BT(new BlockMatrix());
            singleAm1BT->resize(nsrows, nscols);
            if (nscols > 1)
                throw new Exception("Coupling matrix should have one block column");

            if (!BT->block(i, j)->isZero())
            {
                shp<DenseMatrix> BTcollapsed = blockMatrixToDenseMatrix(convert<BlockMatrix>(BT->block(i, j)));

                unsigned int currows = BTcollapsed->nRows();
                unsigned int curcols = BTcollapsed->nCols();

                if (currows != v_dim + p_dim)
                {
                    shp<DENSEMATRIX> inMatrix(new DENSEMATRIX());
                    unsigned int numBlock;
                    if (currows == v_dim){
                        inMatrix->Reshape(p_dim, curcols);
                        numBlock = 1;
                    }
                    else if (currows == p_dim) {
                        inMatrix->Reshape(v_dim, curcols);
                        numBlock = 0;
                    }
                    shp<DenseMatrix> BTnewBlock(new DenseMatrix());
                    BTnewBlock->setMatrix(inMatrix);

                    convert<BlockMatrix>(BT->block(i, j))->setBlock(numBlock, 0,BTnewBlock);

                    currows = v_dim + p_dim;

                    BTcollapsed = blockMatrixToDenseMatrix(convert<BlockMatrix>(BT->block(i, j)));
                }

                std::vector<unsigned int> nrows;
                for (unsigned int irow = 0; irow < nsrows; irow++)
                    nrows.push_back(convert<BlockMatrix>(BT->block(i, j))->block(irow, 0)->nRows());

                for (unsigned int irow = 0; irow < nsrows; irow++)
                {
                    shp<DENSEMATRIX> newMat(new DENSEMATRIX(nrows[irow], curcols));
                    singleAm1BT->setBlock(irow, 0, wrap(newMat));
                }

                for (unsigned int jj = 0; jj < curcols; jj++)
                {
                    // select column vector
                    shp<DENSEVECTOR> curVector(new DENSEVECTOR(currows));
                    shp<DENSEVECTOR> resVector(new DENSEVECTOR(currows));

                    for (unsigned int ii = 0; ii < currows; ii++)
                        (*curVector)(ii) = (*spcast<DENSEMATRIX>(BTcollapsed->data()))(ii, jj);

                    M_solversAs[i]->SetVectors(*resVector, *curVector);
                    M_solversAs[i]->Solve();

                    /*std::cout<<spcast<DENSEMATRIX>(M_collapsedAs[i]->data())->NormInf()<<std::endl<<std::flush;
                    std::cout<<curVector->NormInf()<<std::endl<<std::flush;
                    std::cout<<resVector->NormInf()<<std::endl<<std::flush;*/

                    shp<DENSEVECTOR> a(new DENSEVECTOR(currows));
                    spcast<DENSEMATRIX>(collapsedA->data())->Multiply(false, *resVector, *a);
                    a->Scale(-1.0);
                    a->operator+=(*curVector);

                    shp<DenseVector> curV(new DenseVector());
                    curV->setVector(curVector);
                    shp<DenseVector> resV(new DenseVector());
                    resV->setVector(resVector);

                    shp<DenseVector> rrrVector;
                    rrrVector = convert<DenseVector>(M_collapsedAs[i]->multiplyByVector(resV));
                    curV->multiplyByScalar(-1);
                    rrrVector->add(curV);

                    unsigned int offset = 0;
                    // fill singleAm1BT
                    for (unsigned int irow = 0; irow < nsrows; irow++)
                    {
                        for (unsigned int ii = 0; ii < nrows[irow]; ii++)
                        {
                            (*spcast<DENSEMATRIX>(singleAm1BT->block(irow, 0)->data()))(ii, jj) = (*resVector)(
                                    ii + offset);
                        }
                        offset += nrows[irow];
                    }

                }
            }
            Am1BT->setBlock(i, j, singleAm1BT);
        }

    }

    shp<BlockMatrix> schurComplement(new BlockMatrix());
    schurComplement = convert<BlockMatrix>(B->multiplyByMatrix(Am1BT));
    schurComplement->multiplyByScalar(-1.0);
    schurComplement->add(C);

    M_schurComplementColl = blockMatrixToDenseMatrix(schurComplement);
    M_schurSolver.SetMatrix(*spcast<DENSEMATRIX>(M_schurComplementColl->data()));
}

void
LinearSystemSolver::
solveDense(const BM& matrix, const BV& rhs, BV& sol)
{
    Chrono chrono;
    chrono.start();
    std::string msg = "[LinearSystemSolver] solving dense linear system ...";
    printlog(YELLOW, msg, this->M_data.getVerbose());

    bool monolithicSolve = M_data("rb/online/monolithicSolve", false);

     if (monolithicSolve)  // CURRENTLY NOT WORKING!
     {
         Epetra_SerialDenseSolver monolithicSolver;

         shp<DenseVector> rhsCollapsed = blockVectorToDenseVector(convert<BlockVector>(rhs));
         // rhsCollapsed->dump("rhs");

         msg = " size of linear system = ";
         msg += std::to_string(rhsCollapsed->nRows());
         msg += " ...\n";
         printlog(YELLOW, msg, this->M_data.getVerbose());

         convert<BlockMatrix>(matrix)->printPattern();
         shp<DenseMatrix> matrixCollapsed = blockMatrixToDenseMatrix(convert<BlockMatrix>(matrix));
         DENSEVECTOR solCollapsedDense(rhsCollapsed->nRows());
         // matrixCollapsed->dump("matrix");
         // std::cout << "matrixCollapsed " << spcast<DENSEMATRIX>(matrixCollapsed->data())->NormInf() << std::endl << std::flush;
         // std::cout << "rhsCollapsed " << spcast<DENSEVECTOR>(rhsCollapsed->data())->NormInf() << std::endl << std::flush;
         // std::cout << "nRows " << rhsCollapsed->nRows() << std::endl << std::flush;
         // std::cout << spcast<DENSEMATRIX>(matrixCollapsed->data())->NormInf() << std::endl << std::flush;

         DENSEMATRIX mm = *spcast<DENSEMATRIX>(matrixCollapsed->data());
         DENSEVECTOR vv = *spcast<DENSEVECTOR>(rhsCollapsed->data());
         // std::cout << "matrixCollapsed " << mm.NormInf() << std::endl << std::flush;
         // std::cout << "rhsCollapsed " << vv.NormInf() << std::endl << std::flush;

         monolithicSolver.SetMatrix(mm);
         // std::cout << "code = " << monolithicSolver.Factor() << std::endl << std::flush;

         monolithicSolver.SetVectors(solCollapsedDense, vv);
         std::cout << "Solver Status: " << monolithicSolver.Solve() << std::endl << std::flush;
         // std::cout << spcast<DENSEMATRIX>(matrixCollapsed->data())->NormInf() << std::endl << std::flush;
         // std::cout << "solCollapsed " << solCollapsedDense.NormInf() << std::endl << std::flush;

         exit(1);

         // we need to recast the solCollapsed into a block vector with the same structure
         // as the matrix
         // std::cout << "sol norm = " << sol->norm2() << std::endl << std::flush;
         // shp<BlockVector> solBV = convert<BlockVector>(sol);
         // std::cout << "sol norm = " << solCollapsedDense->Norm2() << std::endl << std::flush;
         // convertVectorType(convert<BlockMatrix>(matrix), wrap(solCollapsedDense), solBV);
         // std::cout << "sol norm = " << solBV->norm2() << std::endl << std::flush;
         // std::cout << "sol norm = " << sol->norm2() << std::endl << std::flush;
         // wrap(solCollapsedDense)->dump("sol");
         // shp<aVector> aux = matrix->multiplyByVector(sol);
         // aux->multiplyByScalar(-1);
         // aux->add(rhs);
         // std::cout << "residual = " << aux->norm2() << std::endl << std::flush;
     }
     else
     {
         bool recomputeSchur = M_data("rb/online/recomputeSchur", false);
         // count primal blocks
         unsigned int nPrimal = 1;
         unsigned int nBlocks = matrix->nRows();

         while (matrix->block(0, nPrimal)->isZero())
             nPrimal++;

         shp<BlockMatrix> A = convert<BlockMatrix>(matrix)->getSubmatrix(0, nPrimal-1, 0, nPrimal-1);
         shp<BlockMatrix> BT = convert<BlockMatrix>(matrix)->getSubmatrix(0, nPrimal-1, nPrimal, nBlocks-1);
         shp<BlockMatrix> B = convert<BlockMatrix>(matrix)->getSubmatrix(nPrimal, nBlocks-1, 0, nPrimal-1);
         shp<BlockMatrix> C = convert<BlockMatrix>(matrix)->getSubmatrix(nPrimal, nBlocks-1, nPrimal, nBlocks-1);

         shp<BlockVector> rhsU = convert<BlockVector>(rhs)->getSubvector(0,nPrimal-1);
         shp<BlockVector> rhsL = convert<BlockVector>(rhs)->getSubvector(nPrimal, nBlocks-1);

         if (recomputeSchur || M_numSolves == 0)
             computeSchurComplementDense(matrix);

         // compute Am1 ru
         shp<BlockVector> Am1ru(new BlockVector());
         Am1ru->resize(nPrimal);

         for (unsigned int i = 0; i < nPrimal; i++)
         {
             shp<DenseVector> collapsedRui = blockVectorToDenseVector(convert<BlockVector>(rhsU->block(i)));
             shp<DENSEVECTOR> solTmp(new DENSEVECTOR(collapsedRui->nRows()));

             M_solversAs[i]->SetVectors(*solTmp, *spcast<DENSEVECTOR>(collapsedRui->data()));
             M_solversAs[i]->Solve();

             unsigned int nrows = rhsU->block(i)->nRows();
             convert<BlockVector>(Am1ru)->setBlock(i,shp<BlockVector>(new BlockVector(nrows)));

             unsigned int offset = 0;
             for (unsigned int ii = 0; ii < nrows; ii++)
             {
                 shp<DENSEVECTOR> dVector(new DENSEVECTOR(rhsU->block(i)->block(ii)->nRows()));

                 for (unsigned int iii = 0; iii < rhsU->block(i)->block(ii)->nRows(); iii++)
                     (*dVector)(iii) = (*solTmp)(iii + offset);

                 offset += rhsU->block(i)->block(ii)->nRows();
                 convert<BlockVector>(Am1ru->block(i))->setBlock(ii, wrap(dVector));
             }
         }

         shp<BlockVector> rhsSchur(new BlockVector());
         rhsSchur->deepCopy(rhsL);
         shp<BlockVector> aux = convert<BlockVector>(B->multiplyByVector(Am1ru));
         aux->multiplyByScalar(-1.0);
         rhsSchur->add(aux);
         shp<DenseVector> rhsSchurCollapsed = blockVectorToDenseVector(rhsSchur);
         shp<DENSEVECTOR> dSolLCollapsed(new DENSEVECTOR(rhsSchurCollapsed->nRows()));
         M_schurSolver.SetVectors(*dSolLCollapsed, *spcast<DENSEVECTOR>(rhsSchurCollapsed->data()));
         M_schurSolver.Solve();

         shp<BlockVector> solL(new BlockVector());
         convertVectorType(B, wrap(dSolLCollapsed), solL);

         shp<BlockVector> rhsA(new BlockVector());
         rhsA->deepCopy(rhsU);
         aux = convert<BlockVector>(BT->multiplyByVector(solL));
         aux->multiplyByScalar(-1.0);
         rhsA->add(aux);

         // solve for u
         shp<BlockVector> solU(new BlockVector());
         solU->resize(nPrimal);

         for (unsigned int i = 0; i < nPrimal; i++)
         {
             shp<DenseVector> collapsedRhsAi = blockVectorToDenseVector(convert<BlockVector>(rhsA->block(i)));
             shp<DENSEVECTOR> solTmp(new DENSEVECTOR(collapsedRhsAi->nRows()));

             M_solversAs[i]->SetVectors(*solTmp, *spcast<DENSEVECTOR>(collapsedRhsAi->data()));
             M_solversAs[i]->Solve();

             shp<DenseMatrix> collapsedA = blockMatrixToDenseMatrix(convert<BlockMatrix>(A->block(i,i)));
             shp<DENSEVECTOR> a(new DENSEVECTOR());
             spcast<DENSEMATRIX>(collapsedA->data())->Multiply(false, *solTmp, *a);
             a->Scale(-1.0);
             a->operator+=(*spcast<DENSEVECTOR>(collapsedRhsAi->data()));

             unsigned int nrows = rhsU->block(i)->nRows();

             shp<BlockVector> iVector(new BlockVector(nrows));

             unsigned int offset = 0;
             for (unsigned int ii = 0; ii < nrows; ii++)
             {
                 shp<DENSEVECTOR> cVector(new DENSEVECTOR(rhsU->block(i)->block(ii)->nRows()));
                 // solU.block(i).block(ii).data().reset(new DENSEVECTOR(rhsU.block(i).block(ii).getNumRows()));

                 for (unsigned int iii = 0; iii < rhsU->block(i)->block(ii)->nRows(); iii++)
                    (*cVector)(iii) = (*solTmp)(iii + offset);

                 offset += rhsU->block(i)->block(ii)->nRows();
                 iVector->setBlock(ii, wrap(cVector));
             }
             solU->setBlock(i, iVector);
         }

         shp<BlockVector> solBV = convert<BlockVector>(sol);
         solBV->resize(rhs->nRows());

         for (unsigned int i = 0; i < solU->nRows(); i++)
             solBV->setBlock(i, solU->block(i));

         unsigned int offset = solU->nRows();
         for (unsigned int i = 0; i < solL->nRows(); i++)
             solBV->setBlock(i + offset, solL->block(i));
     }
     M_numSolves++;

     msg = "done, in ";
     msg += std::to_string(chrono.diff());
     msg += " seconds\n";
     printlog(YELLOW, msg, this->M_data.getVerbose());
 }

void
LinearSystemSolver::
convertVectorType(const shp<BlockMatrix>& matrix,
               const shp<DenseVector>& vector,
               shp<BlockVector>& targetVector)
{
    targetVector->resize(matrix->nRows());

    unsigned int offset = 0;
    for (unsigned int i = 0; i < matrix->nRows(); i++)
    {
     shp<BlockVector> iVector(new BlockVector(matrix->block(i,0)->nRows()));
     for (unsigned int ii = 0; ii < matrix->block(i,0)->nRows(); ii++)
     {
         shp<DENSEVECTOR> inVector(new DENSEVECTOR(matrix->block(i,0)->block(ii,0)->nRows()));
         for (unsigned int iii = 0; iii < matrix->block(i,0)->block(ii,0)->nRows(); iii++)
         {
             (*inVector)(iii) = (*vector)(iii + offset);
         }
         offset += matrix->block(i,0)->block(ii,0)->nRows();
         iVector->setBlock(ii,wrap(inVector));
     }
     targetVector->setBlock(i, iVector);
    }
}

void
LinearSystemSolver::
setPressureMass(const BM &mass)
{
    M_Mp.reset(new BlockMatrix());
    M_Mp->deepCopy(mass);
}

}

