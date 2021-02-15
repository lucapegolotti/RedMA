#include "SaddlePointPreconditioner.hpp"

namespace RedMA
{

SaddlePointPreconditioner::
SaddlePointPreconditioner(const DataContainer& data, const BM& matrix) :
  M_data(data),
  M_matrix(matrix),
  M_thresholdSizeExactSolve(-1)
{
    setup(matrix, true);
}

void
SaddlePointPreconditioner::
setup(const BM& matrix, bool doComputeSchurComplement)
{

    Chrono chrono;
    chrono.start();
    printlog(MAGENTA, "[SaddlePointPreconditioner] starting setup ...\n", M_data.getVerbose());
    // split matrices
    unsigned int nBlocks = matrix->nRows();
    // read options
    setSolverOptions();

    M_maps.reset(new BlockMaps(matrix));
    M_rangeMaps = M_maps->M_rangeMaps;
    M_domainMaps = M_maps->M_domainMaps;
    M_matrixCollapsed = M_maps->M_collapsedMatrix;

    // preconditioner used to approximate Am1 in iterations
    M_innerPrecType = M_data("preconditioner/inner", "SIMPLE");
    M_thresholdSizeExactSolve = M_data("preconditioner/thresholdsize", -1);
    // solution method to approximate BAm1BT in schur
    M_approxSchurType = M_data("preconditioner/approxshur", "SIMPLE");

    unsigned int nPrimal = 1;
    // if nBlocks == 2, then we are imposing weak boundary conditions and
    // we know that nprimal = 1
    if (nBlocks > 2)
    {
        while (nPrimal < nBlocks && matrix->block(0,nPrimal)->isZero())
            nPrimal++;

        /*if (nPrimal == nBlocks)
            throw new Exception("The system does not have a saddle point structure");*/
    }
    unsigned int nDual = nBlocks - nPrimal;
    M_nPrimalBlocks = nPrimal;
    M_nDualBlocks = nDual;

    BM A = matrix->getSubmatrix(0, nPrimal-1, 0, nPrimal-1);
    BM BT(new BlockMatrix(0,0));
    BM B(new BlockMatrix(0,0));
    BM C(new BlockMatrix(0,0));

    M_primalMap.reset(new MAPEPETRA());
    for (unsigned int i = 0; i < 2 * nPrimal; i++)
        *M_primalMap += *M_maps->M_rangeMaps[i];
    printlog(YELLOW,"[SaddlePointPreconditioner] primal map size = " +
                    std::to_string((nPrimal == 0) ? 0 : M_primalMap->mapSize()) +
                    "\n", M_data.getVerbose());

    M_dualMap.reset(new MAPEPETRA());

    if (nDual > 0) {
        for (unsigned int i = 0; i < nDual; i++)
            *M_dualMap += *M_maps->M_rangeMaps[i + nPrimal * 2];
        printlog(YELLOW,"[SaddlePointPreconditioner] dual map size = " +
                        std::to_string(M_dualMap->mapSize()) +
                        "\n", M_data.getVerbose());

        BT = matrix->getSubmatrix(0, nPrimal - 1, nPrimal, nBlocks - 1);
        B = matrix->getSubmatrix(nPrimal, nBlocks - 1, 0, nPrimal - 1);
        C = matrix->getSubmatrix(nPrimal, nBlocks - 1, nPrimal, nBlocks - 1);
    }
    else
        printlog(YELLOW,"[SaddlePointPreconditioner] dual map size = " +
                        std::to_string(0) + ". No coupling occurs!\n", M_data.getVerbose());

    // BlockMaps bmaps(*BT);
    // M_primalMap = bmaps.getMonolithicRangeMapEpetra();
    // M_dualMap = bmaps.getMonolithicDomainMapEpetra();
    //
    // BlockMaps allmaps(*matrix);
    // M_rangeMaps = allmaps.getRangeMapsEpetra();
    // M_domainMaps = allmaps.getDomainMapsEpetra();
    // M_monolithicMap = allmaps.getMonolithicRangeMapEpetra();
    // M_setupTime = chrono.diff();
    //
    // std::string msg = "Monolithic map size = ";
    // msg += std::to_string(M_monolithicMap->mapSize());
    // msg += "\n";
    // printlog(GREEN, msg, M_data.getVerbose());
    //
    // // this is to be optimized
    // M_matrixCollapsed->shallowCopy(collapseBlocks(matrix, allmaps));

    findSmallBlocks(A);

    allocateInnerPreconditioners(A);
    allocateApproximatedInverses(A);

    if ((!std::strcmp(M_innerPrecType.c_str(), "exact") ||
         !std::strcmp(M_approxSchurType.c_str(), "exact")))
        allocateInverseSolvers(A);

    if ((nDual > 0) && (doComputeSchurComplement))
        computeSchurComplement(A, BT, B, C);

    M_setupTime = chrono.diff();
    std::string msg = "done, in ";
    msg += std::to_string(M_setupTime);
    msg += " seconds\n";
    printlog(MAGENTA, msg, M_data.getVerbose());
}

void
SaddlePointPreconditioner::
allocateApproximatedInverses(const BM& primalMatrix)
{
    unsigned int nrows = primalMatrix->nRows();
    M_approximatedInverses.resize(nrows);
    for (unsigned int i = 0; i < nrows; i++)
    {
        if (M_isSmallBlock[i])
        {
            shp<SparseMatrix> Blockcoll = blockMatrixToSparseMatrix(spcast<BlockMatrix>(primalMatrix->block(i,i)));// = collapseBlocks(*Scoll1, BlockMaps(*Scoll1));

            // create invertible approximated matrix
            M_approximatedInverses[i].reset(new ApproxInv);

            shp<Teuchos::ParameterList> smallMatricesOptions;
            smallMatricesOptions.reset(new
                Teuchos::ParameterList(M_solversOptionsInner->sublist("SmallMatricesInverses")));

            // spcast<MATRIXEPETRA>(Scoll->data())->spy("schur");
            M_approximatedInverses[i]->SetRowMatrix(spcast<MATRIXEPETRA>(Blockcoll->data())->matrixPtr());
            M_approximatedInverses[i]->SetParameterList(*smallMatricesOptions);
            M_approximatedInverses[i]->Compute();
        }
    }
}

void
SaddlePointPreconditioner::
findSmallBlocks(const BM& primalMatrix)
{
    unsigned int nrows = primalMatrix->nRows();
    M_isSmallBlock.resize(nrows);
    for (unsigned int i = 0; i < nrows; i++)
    {
        M_isSmallBlock[i] = false;

        unsigned int size1 = spcast<MATRIXEPETRA>(convert<BlockMatrix>(primalMatrix->block(i,i))->block(0,0)->data())->rangeMapPtr()->mapSize();
        unsigned int size2 = spcast<MATRIXEPETRA>(convert<BlockMatrix>(primalMatrix->block(i,i))->block(1,0)->data())->rangeMapPtr()->mapSize();

        if (size1 + size2 < M_thresholdSizeExactSolve)
            M_isSmallBlock[i] = true;
    }
}

void
SaddlePointPreconditioner::
computeSchurComplement(const BM& A, const BM& BT, const BM& B, const BM& C)
{
    computeAm1BT(A, BT);

    M_S = convert<BlockMatrix>(B->multiplyByMatrix(M_Am1BT));
    M_S->multiplyByScalar(-1.0);
    M_S->add(C);

    shp<SparseMatrix> Scoll = blockMatrixToSparseMatrix(M_S);// = collapseBlocks(*Scoll1, BlockMaps(*Scoll1));

    // create invertible approximated matrix
    M_approximatedSchurInverse.reset(new ApproxInv);

    shp<Teuchos::ParameterList> globalSchurOptions;
    globalSchurOptions.reset(new
        Teuchos::ParameterList(M_solversOptionsInner->sublist("GlobalSchurOperator")));

    // spcast<MATRIXEPETRA>(Scoll->data())->spy("schur");
    M_approximatedSchurInverse->SetRowMatrix(spcast<MATRIXEPETRA>(Scoll->data())->matrixPtr());
    M_approximatedSchurInverse->SetParameterList(*globalSchurOptions);
    M_approximatedSchurInverse->Compute();
}

void
SaddlePointPreconditioner::
allocateInverseSolvers(const BM& primalMatrix)
{
    unsigned int nrows = primalMatrix->nRows();
    M_invOperators.resize(nrows);
    for (unsigned int i = 0; i < nrows; i++)
    {
        if (!M_isSmallBlock[i])
        {
            shp<InvOp> invOper;
            shp<NSOp> oper(new NSOp);

            boost::numeric::ublas::matrix<shp<MATRIXEPETRA::matrix_type>> matrixGrid(2,2);
            matrixGrid(0,0) = spcast<MATRIXEPETRA>(convert<BlockMatrix>(primalMatrix->block(i,i))->block(0,0)->data())->matrixPtr();
            matrixGrid(1,0) = spcast<MATRIXEPETRA>(convert<BlockMatrix>(primalMatrix->block(i,i))->block(1,0)->data())->matrixPtr();
            matrixGrid(0,1) = spcast<MATRIXEPETRA>(convert<BlockMatrix>(primalMatrix->block(i,i))->block(0,1)->data())->matrixPtr();
            if (!convert<BlockMatrix>(primalMatrix->block(i,i))->block(1,1)->isZero())
                matrixGrid(1,1) = spcast<MATRIXEPETRA>(convert<BlockMatrix>(primalMatrix->block(i,i))->block(1,1)->data())->matrixPtr();

            unsigned int size1 = spcast<MATRIXEPETRA>(convert<BlockMatrix>(primalMatrix->block(i,i))->block(0,0)->data())->rangeMapPtr()->mapSize();
            unsigned int size2 = spcast<MATRIXEPETRA>(convert<BlockMatrix>(primalMatrix->block(i,i))->block(1,0)->data())->rangeMapPtr()->mapSize();

            if ((M_thresholdSizeExactSolve == -1) || (size1 + size2 < M_thresholdSizeExactSolve))
            {
                oper->setUp(matrixGrid,
                spcast<MATRIXEPETRA>(convert<BlockMatrix>(primalMatrix->block(i,i))->block(0,0)->data())->rangeMap().commPtr());

                shp<Teuchos::ParameterList> options;

                options.reset(new Teuchos::ParameterList(M_solversOptionsInner->sublist("InnerBlockOperator")));

                double innertol = M_data("preconditioner/innertol", 1e-3);

                std::string solverType(options->get<std::string>("Linear Solver Type"));
                if (!std::strcmp(solverType.c_str(),"AztecOO"))
                    options->sublist(solverType).sublist("options")
                                                .get<double>("tol") = innertol;
                else if (!std::strcmp(solverType.c_str(),"Belos"))
                    options->sublist(solverType).sublist("options")
                                                .get<double>("Convergence Tolerance") = innertol;

                invOper.reset(LifeV::Operators::InvertibleOperatorFactory::instance().createObject(solverType));

                invOper->setParameterList(options->sublist(solverType));
                invOper->setOperator(oper);
                invOper->setPreconditioner(M_innerPreconditioners[i]);
                M_invOperators[i] = invOper;
            }
        }
    }
}

void
SaddlePointPreconditioner::
allocateInnerPreconditioners(const BM& primalMatrix)
{
    unsigned int nrows = primalMatrix->nRows();

    M_innerPreconditioners.resize(nrows);

    for (unsigned int i = 0; i < nrows; i++)
    {
        if (!M_isSmallBlock[i])
        {
            if (convert<BlockMatrix>(primalMatrix->block(i,i))->block(0,0)->type() == SPARSE)
            {
                shp<NSPrec> newPrec;
                newPrec.reset(LifeV::Operators::NSPreconditionerFactory::
                              instance().createObject("SIMPLE"));
                Teuchos::RCP<Teuchos::ParameterList> curList(M_solversOptionsInner);
                unsigned int size = spcast<MATRIXEPETRA>(convert<BlockMatrix>(primalMatrix->block(i,i))->block(0,0)->data())->rangeMapPtr()->mapSize();
                // TODO: this obviously needs to be fixed..
                if (size < 10000)
                {
                    curList->sublist("MomentumOperator").get<std::string>("preconditioner type") = "Ifpack";
                }
                else
                {
                    curList->sublist("MomentumOperator").get<std::string>("preconditioner type") = "ML";
                }

                newPrec->setOptions(*curList);

                if (!convert<BlockMatrix>(primalMatrix->block(i,i))->block(1,1)->isZero())
                {
                    newPrec->setUp(spcast<MATRIXEPETRA>(convert<BlockMatrix>(primalMatrix->block(i,i))->block(0,0)->data()),
                                   spcast<MATRIXEPETRA>(convert<BlockMatrix>(primalMatrix->block(i,i))->block(1,0)->data()),
                                   spcast<MATRIXEPETRA>(convert<BlockMatrix>(primalMatrix->block(i,i))->block(0,1)->data()),
                                   spcast<MATRIXEPETRA>(convert<BlockMatrix>(primalMatrix->block(i,i))->block(1,1)->data()));
                }
                else
                {
                    newPrec->setUp(spcast<MATRIXEPETRA>(convert<BlockMatrix>(primalMatrix->block(i,i))->block(0,0)->data()),
                                   spcast<MATRIXEPETRA>(convert<BlockMatrix>(primalMatrix->block(i,i))->block(1,0)->data()),
                                   spcast<MATRIXEPETRA>(convert<BlockMatrix>(primalMatrix->block(i,i))->block(0,1)->data()));
                }

                std::vector<shp<Epetra_Map>> localRangeMaps(2);
                std::vector<shp<Epetra_Map>> localDomainMaps(2);

                for (unsigned int j = 0; j < 2; j++)
                    localRangeMaps[j].reset(new Epetra_Map(
                    spcast<MATRIXEPETRA>(convert<BlockMatrix>(primalMatrix->block(i,i))->block(j,0)->data())->matrixPtr()->OperatorRangeMap()));

                for (unsigned int j = 0; j < 2; j++)
                    localDomainMaps[j].reset(new Epetra_Map(
                    spcast<MATRIXEPETRA>(convert<BlockMatrix>(primalMatrix->block(i,i))->block(0,j)->data())->matrixPtr()->OperatorDomainMap()));

                shp<LifeV::BlockEpetra_Map> rmblock(new LifeV::BlockEpetra_Map(localRangeMaps));
                shp<LifeV::BlockEpetra_Map> dmblock(new LifeV::BlockEpetra_Map(localDomainMaps));
                newPrec->setRangeMap(rmblock);
                newPrec->setDomainMap(dmblock);

                newPrec->updateApproximatedMomentumOperator();
                newPrec->updateApproximatedSchurComplementOperator();
                M_innerPreconditioners[i] = newPrec;
            }
        }
    }
}

void
SaddlePointPreconditioner::
computeAm1BT(const BM& A, const BM& BT)
{
    M_Am1BT.reset(new BlockMatrix(A->nRows(), BT->nCols()));

    for (unsigned int i = 0; i < A->nRows(); i++)
    {
        for (unsigned int j = 0; j < BT->nCols(); j++)
        {
            BM Abm = spcast<BlockMatrix>(A->block(i,i));
            BM BTbm = spcast<BlockMatrix>(BT->block(i,j));
            M_Am1BT->setBlock(i,j,computeSingleAm1BT(Abm, BTbm,i));
        }
    }
}

shp<aMatrix>
SaddlePointPreconditioner::
computeSingleAm1BT(const BM& A, const BM& BT,
                   const unsigned int& index)
{
    shp<BlockMatrix> retMat;

    if (!BT->isZero())
    {
        retMat.reset(new BlockMatrix(A->nRows(), BT->nCols()));
        Chrono chrono;
        chrono.start();
        printlog(GREEN, "[SaddlePointPreconditioner] single AM1BT ...", M_data.getVerbose());

        MAPEPETRA rangeMapU = static_cast<MATRIXEPETRA*>(BT->block(0,0)->data().get())->rangeMap();
        MAPEPETRA rangeMapP = static_cast<MATRIXEPETRA*>(A->block(1,0)->data().get())->rangeMap();
        MAPEPETRA domainMap = static_cast<MATRIXEPETRA*>(BT->block(0,0)->data().get())->domainMap();

        unsigned int ncols = domainMap.mapSize();

        std::vector<shp<DistributedVector>> ressU(ncols);
        std::vector<shp<DistributedVector>> ressP(ncols);

        bool exactSolve = false;

        if (!std::strcmp(M_approxSchurType.c_str(),"exact") ||
            M_isSmallBlock[index]) exactSolve = true;

        VECTOREPETRA selector(domainMap);
        VECTOREPETRA colU(rangeMapU);
        VECTOREPETRA colP(rangeMapP);

        for (unsigned int i = 0; i < ncols; i++)
        {
            selector.zero();
            colU.zero();
            colP.zero();
            if (domainMap.isOwned(i))
                selector[i] = 1.0;

            colU = (*spcast<MATRIXEPETRA>(BT->block(0,0)->data())) * selector;

            ressU[i].reset(new DistributedVector());
            ressU[i]->setData(shp<VECTOREPETRA>(new VECTOREPETRA(rangeMapU)));
            ressP[i].reset(new DistributedVector());
            ressP[i]->setData(shp<VECTOREPETRA>(new VECTOREPETRA(rangeMapP)));

            if (exactSolve)
            {
                auto monolithicMap = rangeMapU;
                monolithicMap += rangeMapP;

                VECTOREPETRA X(monolithicMap);
                VECTOREPETRA Y(monolithicMap);
                X.zero();
                Y.zero();

                X.subset(colU, rangeMapU, 0, 0);
                CoutRedirecter ct;
                ct.redirect();
                if (M_isSmallBlock[index])
                    M_approximatedInverses[index]->ApplyInverse(X.epetraVector(), Y.epetraVector());
                else
                    M_invOperators[index]->ApplyInverse(X.epetraVector(), Y.epetraVector());
                ct.restore();

                spcast<VECTOREPETRA>(ressU[i]->data())->subset(Y, rangeMapU, 0, 0);
                spcast<VECTOREPETRA>(ressP[i]->data())->subset(Y, rangeMapP, rangeMapU.mapSize(), 0);
            }
            else
            {
                M_innerPreconditioners[index]->ApplyInverse(colU, colP,
                                                            *spcast<VECTOREPETRA>(ressU[i]->data()),
                                                            *spcast<VECTOREPETRA>(ressP[i]->data()));
            }
        }

        shp<SparseMatrix> Am1BTu(new SparseMatrix(ressU));
        shp<SparseMatrix> Am1BTp(new SparseMatrix(ressP));
        retMat->setBlock(0,0,Am1BTu);
        retMat->setBlock(1,0,Am1BTp);

        std::string msg = "done, in ";
        msg += std::to_string(chrono.diff());
        msg += " seconds\n";
        printlog(GREEN, msg, M_data.getVerbose());
    }
    return retMat;
}

void
SaddlePointPreconditioner::
setSolverOptions()
{
    std::string optionsPrec = M_data("preconditioner/options",
                                     "datafiles/solversOptionsFast");
    optionsPrec += ".xml";
    M_solversOptionsInner = Teuchos::getParametersFromXmlFile(optionsPrec);
    shp<Teuchos::ParameterList> monolithicOptions;
    monolithicOptions.reset(
        new Teuchos::ParameterList(M_solversOptionsInner->sublist("InnerBlockOperator")));
    M_pListLinSolver = monolithicOptions;
}

void
SaddlePointPreconditioner::
solveEveryPrimalBlock(const VECTOREPETRA& X, VECTOREPETRA &Y) const
{
    unsigned int offset = 0;
    unsigned int count = 0;
    for (unsigned int i = 0; i < M_nPrimalBlocks; i++)
    {
        MAPEPETRA rangeVelocity = *M_rangeMaps[count];
        MAPEPETRA rangePressure = *M_rangeMaps[count+1];

        unsigned int size1 = rangeVelocity.mapSize();
        unsigned int size2 = rangePressure.mapSize();

        bool exactSolve = false;

        if (!std::strcmp(M_innerPrecType.c_str(),"exact") ||
            M_isSmallBlock[i]) exactSolve = true;

        if (exactSolve)
        {
            auto monolithicMap = rangeVelocity;
            monolithicMap += rangePressure;

            VECTOREPETRA subX(monolithicMap);
            VECTOREPETRA subY(monolithicMap);

            subX.subset(X, monolithicMap, offset, 0);

            RedMA::CoutRedirecter ct;
            ct.redirect();

            if (M_isSmallBlock[i])
                M_approximatedInverses[i]->ApplyInverse(subX.epetraVector(),
                                                        subY.epetraVector());
            else
                M_invOperators[i]->ApplyInverse(subX.epetraVector(),
                                                subY.epetraVector());

            ct.restore();

            Y.subset(subY, monolithicMap, 0, offset);
        }
        else
        {
            VECTOREPETRA subVelocity(rangeVelocity);
            VECTOREPETRA subPressure(rangePressure);
            subVelocity.zero();
            subPressure.zero();

            subVelocity.subset(X, rangeVelocity, offset, 0);
            subPressure.subset(X, rangePressure,
                               offset + rangeVelocity.mapSize(), 0);

            VECTOREPETRA resVelocity(rangeVelocity);
            VECTOREPETRA resPressure(rangePressure);

            M_innerPreconditioners[i]->ApplyInverse(subVelocity, subPressure,
                                                    resVelocity, resPressure);

            Y.subset(resVelocity, rangeVelocity, 0, offset);
            Y.subset(resPressure, rangePressure, 0,
                     offset + rangeVelocity.mapSize());
        }
        offset += rangeVelocity.mapSize() + rangePressure.mapSize();
        count = count + 2;
    }
}

void
SaddlePointPreconditioner::
applyEveryB(const VECTOREPETRA& X, VECTOREPETRA &Y) const
{
    unsigned int offset = 0;

    std::vector<VECTOREPETRA> Xs;
    for (unsigned int i = 0; i < M_nPrimalBlocks * 2; i++)
    {
        VECTOREPETRA newVec(*M_rangeMaps[i]);
        newVec.subset(X, *M_rangeMaps[i], offset, 0);
        Xs.push_back(newVec);
        offset += M_rangeMaps[i]->mapSize();
    }

    offset = 0;
    for (unsigned int i = M_nPrimalBlocks * 2; i < M_nPrimalBlocks * 2 + M_nDualBlocks; i++)
    {
        MAPEPETRA curRange = *M_rangeMaps[i];
        VECTOREPETRA subRes(curRange, LifeV::Unique);
        subRes.zero();
        for (unsigned int j = 0; j < M_nPrimalBlocks * 2; j++)
        {
            if (!M_matrixCollapsed->block(i,j)->isZero())
                subRes += (*spcast<MATRIXEPETRA>(M_matrixCollapsed->block(i,j)->data())) * Xs[j];
        }
        Y.subset(subRes, curRange, 0, offset);
        offset += curRange.mapSize();
    }
}

void
SaddlePointPreconditioner::
applyEveryAm1BT(const VECTOREPETRA& X, VECTOREPETRA &Y) const
{
    unsigned int offset = 0;

    std::vector<VECTOREPETRA> Xs;
    for (unsigned int i = M_nPrimalBlocks * 2; i < M_nPrimalBlocks * 2 + M_nDualBlocks; i++)
    {
        VECTOREPETRA newVec(*M_domainMaps[i]);
        newVec.subset(X, *M_domainMaps[i], offset, 0);
        Xs.push_back(newVec);
        offset += M_domainMaps[i]->mapSize();
    }

    offset = 0;
    for (unsigned int i = 0; i < 2 * M_nPrimalBlocks; i++)
    {
        MAPEPETRA curRange = *M_domainMaps[i];
        VECTOREPETRA subRes(curRange, LifeV::Unique);
        subRes.zero();
        // for (unsigned int j = 2 * M_nPrimalBlocks; j < M_nPrimalBlocks * 2 + M_nDualBlocks; j++)
        // {
        //     if (!M_matrixCollapsed.block(i,j).isNull())
        //     {
        //         subRes += (*M_matrixCollapsed.block(i,j).data()) * Xs[j-2*M_nPrimalBlocks];
        //     }
        // }
        for (unsigned int j = 0; j < M_nDualBlocks; j++)
        {
            if (!M_Am1BT->block(i/2,j)->isZero())
                if (!convert<BlockMatrix>(M_Am1BT->block(i/2,j))->block(i%2,0)->isZero())
                    subRes += (*spcast<MATRIXEPETRA>(convert<BlockMatrix>(M_Am1BT->block(i/2,j))->block(i%2,0)->data())) * Xs[j];
        }
        Y.subset(subRes, curRange, 0, offset);
        offset += curRange.mapSize();
    }
}


int
SaddlePointPreconditioner::
ApplyInverse(const super::vector_Type& X, super::vector_Type& Y) const
{
    if (M_nDualBlocks > 0)
    {
        const VECTOREPETRA X_vectorEpetra(X, M_maps->M_monolithicRangeMap, LifeV::Unique);
        VECTOREPETRA Y_vectorEpetra(M_maps->M_monolithicRangeMap, LifeV::Unique);
        VECTOREPETRA X_primal(M_primalMap, LifeV::Unique);
        VECTOREPETRA X_dual(M_dualMap, LifeV::Unique);
        X_primal.subset(X_vectorEpetra, *M_primalMap, 0, 0);
        X_dual.subset(X_vectorEpetra, *M_dualMap, M_primalMap->mapSize(), 0);

        // here we store the result
        VECTOREPETRA Y_primal(M_primalMap, LifeV::Unique);
        VECTOREPETRA Y_dual(M_dualMap, LifeV::Unique);

        VECTOREPETRA Z(M_primalMap, LifeV::Unique);
        solveEveryPrimalBlock(X_primal, Z);

        VECTOREPETRA Bz(M_dualMap, LifeV::Unique);
        applyEveryB(Z, Bz);

        M_approximatedSchurInverse->ApplyInverse((X_dual-Bz).epetraVector(),
                                                 Y_dual.epetraVector());

        VECTOREPETRA Am1BTy(M_primalMap, LifeV::Unique);
        applyEveryAm1BT(Y_dual, Am1BTy);

        // this can be optimized because we already have Am1BTy
        // solveEveryPrimalBlock(BTy, Am1BTy);

        Y_primal = Z - Am1BTy;
        Y_vectorEpetra.subset(Y_primal, *M_primalMap, 0, 0);
        Y_vectorEpetra.subset(Y_dual, *M_dualMap, 0, M_primalMap->mapSize());
        Y = dynamic_cast<Epetra_MultiVector&>(Y_vectorEpetra.epetraVector());
    }
    else
    {
        M_innerPreconditioners[0]->ApplyInverse(X, Y);
    }
    return 0;
}

}
