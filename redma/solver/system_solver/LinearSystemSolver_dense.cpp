#include "LinearSystemSolver.hpp"

namespace RedMA
{

template <>
void
LinearSystemSolver<BlockVector<DenseVector>, BlockMatrix<DenseMatrix>>::
buildPreconditioner(const BlockMatrix<BlockMatrix<DenseMatrix>>& matrix)
{
}

template <>
void
LinearSystemSolver<BlockVector<DenseVector>, BlockMatrix<DenseMatrix>>::
computeSchurComplementDense(const BM& matrix)
{
    // count primal blocks
    unsigned int nPrimal = 1;
    unsigned int nBlocks = matrix.nRows();

    while (matrix.block(0,nPrimal).isNull())
        nPrimal++;
    BM A = matrix.getSubmatrix(0, nPrimal-1, 0, nPrimal-1);
    BM BT = matrix.getSubmatrix(0, nPrimal-1, nPrimal, nBlocks-1);
    BM B = matrix.getSubmatrix(nPrimal, nBlocks-1, 0, nPrimal-1);
    BM C = matrix.getSubmatrix(nPrimal, nBlocks-1, nPrimal, nBlocks-1);

    M_collapsedAs.resize(nPrimal);
    M_solversAs.resize(nPrimal);

    unsigned int primalIndex = 0;
    for (unsigned int i = 0; i < nPrimal; i++)
    {
        // to do: try to compute factorization and see if it changes anything
        M_solversAs[i].reset(new Epetra_SerialDenseSolver());
        M_collapsedAs[i] = A.block(i,i).collapse();
        M_solversAs[i]->SetMatrix(*M_collapsedAs[i].data());
        M_solversAs[i]->Factor();
    }

    // compute Am1BT
    BM Am1BT;
    Am1BT.resize(BT.nRows(), BT.nCols());
    for (unsigned int i = 0; i < BT.nRows(); i++)
    {
        for (unsigned int j = 0; j < BT.nCols(); j++)
        {
            unsigned int nsrows = BT.block(i,j).nRows();
            unsigned int nscols = BT.block(i,j).nCols();
            BlockMatrix<DenseMatrix> singleAm1BT;
            singleAm1BT.resize(nsrows, nscols);
            if (nscols > 1)
                throw new Exception("Coupling matrix should have one block column");

            std::vector<unsigned int> nrows;
            for (unsigned int irow = 0; irow < nsrows; irow++)
                nrows.push_back(BT.block(i,j).block(irow,0).getNumRows());

            if (!BT.block(i,j).isNull())
            {
                DenseMatrix BTcollapsed = BT.block(i,j).collapse();

                unsigned int currows = BTcollapsed.getNumRows();
                unsigned int curcols = BTcollapsed.getNumCols();

                for (unsigned int irow = 0; irow < nsrows; irow++)
                    singleAm1BT.block(irow,0).data().reset(new DENSEMATRIX(nrows[irow], curcols));

                for (unsigned int jj = 0; jj < curcols; jj++)
                {
                    // select colum vector
                    SHP(DENSEVECTOR) curVector(new DENSEVECTOR(currows));
                    SHP(DENSEVECTOR) resVector(new DENSEVECTOR(currows));

                    for (unsigned int ii = 0 ; ii < currows; ii++)
                        (*curVector)(ii) = (*BTcollapsed.data())(ii,jj);

                    // collapsedAs[i] = A.block(i,i).collapse();
                    // M_solversAs[i]->SetMatrix(*collapsedAs[i].data());
                    M_solversAs[i]->SetVectors(*resVector, *curVector);
                    M_solversAs[i]->Solve();
                    DENSEVECTOR a;
                    A.block(i,i).collapse().data()->Multiply(false, *resVector, a);
                    a.Scale(-1);
                    a += *curVector;
                    DenseVector curV;
                    curV.data() = curVector;
                    DenseVector resV;
                    resV.data() = resVector;

                    DenseVector rrrVector;
                    rrrVector = A.block(i,i).collapse() * resV;
                    rrrVector -= curV;

                    unsigned int offset = 0;
                    // fill singleAm1BT
                    for (unsigned int irow = 0; irow < nsrows; irow++)
                    {
                        for (unsigned int ii = 0; ii < nrows[irow]; ii++)
                        {
                            (*singleAm1BT.block(irow,0).data())(ii,jj) = (*resVector)(ii + offset);
                        }
                        offset += nrows[irow];
                    }
                }
            }
            Am1BT.block(i,j) = singleAm1BT;
        }
    }
    Am1BT.finalize();

    BM schurComplement;
    schurComplement.softCopy(B * Am1BT);
    schurComplement *= (-1.0);
    schurComplement += C;

    M_schurComplementColl = schurComplement.collapse().block(0,0);
    M_schurSolver.SetMatrix(*M_schurComplementColl.data());
}

template <>
void
LinearSystemSolver<BlockVector<DenseVector>, BlockMatrix<DenseMatrix>>::
solve(const BlockMatrix<BlockMatrix<DenseMatrix>>& matrix,
      const BlockVector<BlockVector<DenseVector>>& rhs,
      BlockVector<BlockVector<DenseVector>>& sol)
{
    Chrono chrono;
    chrono.start();

    std::string msg = "[LinearSystemSolver] solving dense linear system ...";
    printlog(YELLOW, msg, this->M_data.getVerbose());

    bool monolithicSolve = M_data("rb/online/monolithicSolve", false);

    if (monolithicSolve)
    {
        Epetra_SerialDenseSolver monolithicSolver;

        DenseVector rhsCollapsed = rhs.collapse().block(0);

        msg = " size of linear system = ";
        msg += std::to_string(rhsCollapsed.getNumRows());
        msg += " ...";
        printlog(YELLOW, msg, this->M_data.getVerbose());

        DenseMatrix matrixCollapsed = matrix.collapse().block(0,0);
        DenseVector solCollapsed;
        solCollapsed.data().reset(new DENSEVECTOR(rhsCollapsed.getNumRows()));

        monolithicSolver.SetMatrix(*matrixCollapsed.data());
        monolithicSolver.SetVectors(*solCollapsed.data(), *rhsCollapsed.data());
        monolithicSolver.Solve();

        matrix.convertVectorType(solCollapsed, sol);
    }
    else
    {
        bool recomputeSchur = M_data("rb/online/recomputeSchur", false);
        // count primal blocks
        unsigned int nPrimal = 1;
        unsigned int nBlocks = matrix.nRows();

        while (matrix.block(0,nPrimal).isNull())
            nPrimal++;

        BM A = matrix.getSubmatrix(0, nPrimal-1, 0, nPrimal-1);
        BM BT = matrix.getSubmatrix(0, nPrimal-1, nPrimal, nBlocks-1);
        BM B = matrix.getSubmatrix(nPrimal, nBlocks-1, 0, nPrimal-1);
        BM C = matrix.getSubmatrix(nPrimal, nBlocks-1, nPrimal, nBlocks-1);

        std::cout << "rhs norm = " << rhs.norm2() << std::endl << std::flush;
        BV rhsU = rhs.getSubvector(0,nPrimal-1);
        std::cout << "rhsU norm = " << rhsU.norm2() << std::endl;
        BV rhsL = rhs.getSubvector(nPrimal, nBlocks-1);
        std::cout << "rhsL norm = " << rhsL.norm2() << std::endl;

        if (recomputeSchur || M_numSolves == 0)
            computeSchurComplementDense(matrix);

        // compute Am1 ru
        BV Am1ru;
        Am1ru.resize(nPrimal);

        for (unsigned int i = 0; i < nPrimal; i++)
        {
            DenseVector collapsedRui = rhsU.block(i).collapse();
            SHP(DENSEVECTOR) sol(new DENSEVECTOR(collapsedRui.getNumRows()));
            M_solversAs[i]->SetVectors(*sol, *collapsedRui.data());
            M_solversAs[i]->Solve();

            unsigned int nrows = rhsU.block(i).nRows();
            Am1ru.block(i).resize(nrows);

            unsigned int offset = 0;
            for (unsigned int ii = 0; ii < nrows; ii++)
            {
                Am1ru.block(i).block(ii).data().reset(new DENSEVECTOR(rhsU.block(i).block(ii).getNumRows()));

                for (unsigned int iii = 0; iii < rhsU.block(i).block(ii).getNumRows(); iii++)
                    (*Am1ru.block(i).block(ii).data())(iii) = (*sol)(iii + offset);

                offset += rhsU.block(i).block(ii).getNumRows();
            }
        }

        BV rhsSchur;
        rhsSchur.hardCopy(rhsL);
        rhsSchur -= B * Am1ru;
        DenseVector rhsSchurCollapsed = rhsSchur.collapse().block(0);
        DenseVector solLCollapsed;
        solLCollapsed.data().reset(new DENSEVECTOR(rhsSchurCollapsed.getNumRows()));
        M_schurSolver.SetVectors(*solLCollapsed.data(), *rhsSchurCollapsed.data());
        M_schurSolver.Solve();
        std::cout << "solLCollapsed norm = " << solLCollapsed.norm2() << std::endl << std::flush;
        BV solL;
        B.convertVectorType(solLCollapsed, solL);

        BV rhsA;
        rhsA.hardCopy(rhsU);
        std::cout << "norm rhsa1 = " << rhsA.norm2() << std::endl << std::flush;
        rhsA -= BT * solL;
        std::cout << "norm rhsa2 = " << rhsA.norm2() << std::endl << std::flush;
        std::cout << "solL norm = " << solL.norm2() << std::endl << std::flush;
        // solve for u
        BV solU;
        solU.resize(nPrimal);

        for (unsigned int i = 0; i < nPrimal; i++)
        {
            DenseVector collapsedRhsAi = rhsA.block(i).collapse();
            SHP(DENSEVECTOR) sol(new DENSEVECTOR(collapsedRhsAi.getNumRows()));

            std::cout << "i = " << i << std::endl << std::flush;
            // collapsedAs[i] = A.block(i,i).collapse();
            // M_solversAs[i]->SetMatrix(*collapsedAs[i].data());
            M_solversAs[i]->SetVectors(*sol, *collapsedRhsAi.data());
            std::cout << "collapsedRhsAi = " << collapsedRhsAi.norm2() << std::endl << std::flush;
            std::cout << "sol = " << sol->Norm2() << std::endl << std::flush;
            M_solversAs[i]->Solve();
            DENSEVECTOR a;

            A.block(i,i).collapse().data()->Multiply(false, *sol, a);
            a.Scale(-1);
            a += *collapsedRhsAi.data();

            std::cout << "residual = " << a.Norm2() << std::endl << std::flush;

            unsigned int nrows = rhsU.block(i).nRows();
            solU.block(i).resize(nrows);

            unsigned int offset = 0;
            for (unsigned int ii = 0; ii < nrows; ii++)
            {
                solU.block(i).block(ii).data().reset(new DENSEVECTOR(rhsU.block(i).block(ii).getNumRows()));

                for (unsigned int iii = 0; iii < rhsU.block(i).block(ii).getNumRows(); iii++)
                    (*solU.block(i).block(ii).data())(iii) = (*sol)(iii + offset);

                offset += rhsU.block(i).block(ii).getNumRows();
            }
        }

        std::cout << "solU norm = " << solU.norm2() << std::endl << std::flush;

        // soft copy sol u in solution vector
        sol.resize(rhs.nRows());
        for (unsigned int i = 0; i < solU.nRows(); i++)
        {
            sol.block(i).resize(solU.block(i).nRows());
            for (unsigned int ii = 0; ii < solU.block(i).nRows(); ii++)
            {
                sol.block(i).block(ii).softCopy(solU.block(i).block(ii));
            }
        }

        unsigned int offset = solU.nRows();
        for (unsigned int i = 0; i < solL.nRows(); i++)
        {
            sol.block(i + offset).resize(solL.block(i).nRows());
            for (unsigned int ii = 0; ii < solL.block(i).nRows(); ii++)
            {
                sol.block(i + offset).block(ii).softCopy(solL.block(i).block(ii));
            }
        }
    }
    M_numSolves++;

    msg = "done, in ";
    msg += std::to_string(chrono.diff());
    msg += " seconds\n";
    printlog(YELLOW, msg, this->M_data.getVerbose());
    // exit(1);
}

template <>
void
LinearSystemSolver<Double, Double>::
solve(const BlockMatrix<Double>& matrix,
      const BlockVector<Double>& rhs,
      BlockVector<Double>& sol)
{
    if (matrix.nRows() > 1 || matrix.nCols() > 1 || rhs.nRows() > 1)
        throw new Exception("solve with blocks of doubles is implemented only in"
                            " dimension 1");

    sol.block(0).data() = rhs.block(0).data() / matrix.block(0,0).data();

    M_numSolves++;
}

}
