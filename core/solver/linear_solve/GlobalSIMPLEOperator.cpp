#include <GlobalSIMPLEOperator.hpp>

#include <lifev/core/linear_algebra/IfpackPreconditioner.hpp>
#include <lifev/core/linear_algebra/MLPreconditioner.hpp>
#include <lifev/core/linear_algebra/TwoLevelPreconditioner.hpp>
#include <lifev/core/linear_algebra/AztecooOperatorAlgebra.hpp>
#include <lifev/core/linear_algebra/BelosOperatorAlgebra.hpp>
#include <lifev/core/linear_algebra/ApproximatedInvertibleRowMatrix.hpp>

#include <lifev/core/algorithm/PreconditionerML.hpp>

namespace LifeV
{
namespace Operators
{

GlobalSIMPLEOperator::
GlobalSIMPLEOperator(std::string singleOperatorType, bool exactSolve,
                     int numIterationsExactSolve) :
  M_label("GlobalSIMPLEOperator"),
  M_useTranspose(false),
  M_singleOperatorType(singleOperatorType),
  M_solvePrimalBlocksExactly(exactSolve),
  M_numIterationExactSolve(numIterationsExactSolve),
  M_countIterations(0)
{

}

GlobalSIMPLEOperator::~GlobalSIMPLEOperator()
{

}

void
GlobalSIMPLEOperator::
setUp(RedMA::BlockMatrix matrix,
      const commPtr_Type & comm)
{
    LifeChrono chrono;
    chrono.start();
    RedMA::printlog(RedMA::MAGENTA, "[GlobalSIMPLEOperator] starting setup ...\n",
                    M_verbose);
    M_matrix = matrix;

    operatorPtrContainer_Type blockOper = matrix.getGrid();

    M_nBlockRows = blockOper.size1();
    M_nBlockCols = blockOper.size2();
    M_comm = comm;

    BlockEpetra_Map::mapPtrContainer_Type rangeBlockMaps(M_nBlockRows);
    BlockEpetra_Map::mapPtrContainer_Type domainBlockMaps(M_nBlockCols);

    M_nPrimalBlocks = 0;
    for (UInt iblock=0; iblock < M_nBlockRows; ++iblock)
    {
        for (UInt jblock=0; jblock < M_nBlockCols; ++jblock)
        {
            if (blockOper(iblock,jblock) != 0 && rangeBlockMaps[iblock]==0)
            {
                rangeBlockMaps[iblock].reset(new
                      Epetra_Map(blockOper(iblock,jblock)->OperatorRangeMap()));
                jblock = M_nBlockCols;
            }
        }
        if (iblock % 2 == 0 && blockOper(iblock,iblock) != 0)
        {
            // compute domain and range maps
            BlockEpetra_Map::mapPtrContainer_Type localRangeBlockMaps(2);
            BlockEpetra_Map::mapPtrContainer_Type localDomainBlockMaps(2);

            localRangeBlockMaps[0].reset(new
                    Epetra_Map(blockOper(iblock,iblock)->OperatorRangeMap()));

            // this could be the case in steady case for example
            if (blockOper(iblock,iblock+1) != 0)
                localRangeBlockMaps[1].reset(new
                        Epetra_Map(blockOper(iblock,iblock+1)->OperatorRangeMap()));
            else
                localRangeBlockMaps[1].reset(new
                        Epetra_Map(blockOper(iblock+1,iblock+1)->OperatorRangeMap()));

            localDomainBlockMaps[0].reset(new
                    Epetra_Map(blockOper(iblock,iblock)->OperatorDomainMap()));

            if (blockOper(iblock+1,iblock) != 0)
                localDomainBlockMaps[1].reset(new
                        Epetra_Map(blockOper(iblock+1,iblock)->OperatorDomainMap()));
            else
                localDomainBlockMaps[1].reset(new
                        Epetra_Map(blockOper(iblock+1,iblock+1)->OperatorDomainMap()));

            M_SingleOperators.push_back(allocateSingleOperator(matrix, iblock,
                                                               localRangeBlockMaps,
                                                               localDomainBlockMaps));

            if (M_solvePrimalBlocksExactly)
            {
                InvOperatorPtr newInvOper = allocateSingleInvOperator(matrix, iblock);
                newInvOper->setPreconditioner(M_SingleOperators[M_nPrimalBlocks]);;
                M_invOperators.push_back(newInvOper);
            }

            M_nPrimalBlocks++;
        }
    }
    for (UInt jblock=0; jblock < M_nBlockCols; ++jblock)
        for (UInt iblock=0; iblock < M_nBlockRows; ++iblock)
        {
            if (blockOper(iblock,jblock) != 0 && domainBlockMaps[jblock]==0)
            {
                domainBlockMaps[jblock].reset(new
                     Epetra_Map(blockOper(iblock,jblock)->OperatorDomainMap()));
                iblock = M_nBlockRows;
            }
        }

    M_domainMap.reset(new BlockEpetra_Map(domainBlockMaps));
    M_rangeMap.reset(new BlockEpetra_Map(rangeBlockMaps));

    M_oper = blockOper;

    fillComplete();
    if (M_nPrimalBlocks > 1)
        computeGlobalSchurComplement();

    std::string msg = "done, in ";
    msg += std::to_string(chrono.diff());
    msg += " seconds\n";
    RedMA::printlog(RedMA::MAGENTA, msg, M_verbose);
}

GlobalSIMPLEOperator::InvOperatorPtr
GlobalSIMPLEOperator::
allocateSingleInvOperator(RedMA::BlockMatrix matrix,
                          UInt iblock)
{
    // here we rely on the structure of the global matrix
    InvOperatorPtr invOper;

    std::shared_ptr<LifeV::Operators::NavierStokesOperator> oper(
                                     new LifeV::Operators::NavierStokesOperator);

    operatorPtrContainer_Type operData(2,2);

    operData(0,0) = matrix.block(iblock,iblock)->matrixPtr();
    operData(1,0) = matrix.block(iblock+1,iblock)->matrixPtr();
    operData(0,1) = matrix.block(iblock,iblock+1)->matrixPtr();

    // then we are using stabilization
    if (matrix.block(iblock+1,iblock+1) != nullptr)
    {
        operData(1,1) = matrix.block(iblock+1,iblock+1)->matrixPtr();
    }

    oper->setUp(operData, M_comm);

    std::shared_ptr<Teuchos::ParameterList> options;
    options.reset(
        new Teuchos::ParameterList(M_solversOptions.sublist("InnerBlockOperator")));

    std::string solverType(options->get<std::string>("Linear Solver Type"));
    invOper.reset(LifeV::Operators::InvertibleOperatorFactory::instance().createObject(solverType));
    invOper->setParameterList(options->sublist(solverType));

    invOper->setOperator(oper);

    return invOper;
}

GlobalSIMPLEOperator::PreconditionerPtr
GlobalSIMPLEOperator::
allocateSingleOperator(RedMA::BlockMatrix matrix, UInt iblock,
                       BlockEpetra_Map::mapPtrContainer_Type localRangeBlockMaps,
                       BlockEpetra_Map::mapPtrContainer_Type localDomainBlockMaps)
{
    // here we rely on the structure of the global matrix
    PreconditionerPtr newPrec;
    if (!std::strcmp(M_singleOperatorType.c_str(), "SIMPLE"))
    {
        newPrec.reset(Operators::NSPreconditionerFactory::
                      instance().createObject("SIMPLE"));
        newPrec->setOptions(M_solversOptions);

        // then we are not considering stabilization
        if (matrix.block(iblock+1,iblock+1) == nullptr)
        {
            newPrec->setUp(matrix.block(iblock,iblock),
                           matrix.block(iblock+1,iblock),
                           matrix.block(iblock,iblock+1));
        }
        // then we are using stabilization
        else
        {
            newPrec->setUp(matrix.block(iblock,iblock),
                           matrix.block(iblock+1,iblock),
                           matrix.block(iblock,iblock+1),
                           matrix.block(iblock+1,iblock+1));
        }

        std::shared_ptr<BlockEpetra_Map> localRangeMap(
                                  new BlockEpetra_Map(localRangeBlockMaps));
        std::shared_ptr<BlockEpetra_Map> localDomainMap(
                                  new BlockEpetra_Map(localDomainBlockMaps));

        newPrec->setRangeMap(localRangeMap);
        newPrec->setDomainMap(localDomainMap);
        newPrec->updateApproximatedMomentumOperator();
        newPrec->updateApproximatedSchurComplementOperator();
    }
    else if (!std::strcmp(M_singleOperatorType.c_str(), "STEADY"))
    {
        newPrec.reset(Operators::NSPreconditionerFactory::
                      instance().createObject("STEADY"));
        newPrec->setOptions(M_solversOptions);
        // newPrec->SetSolveMomentumBelos();

        newPrec->setUp(matrix.block(iblock,iblock),
                       matrix.block(iblock,iblock+1),
                       matrix.block(iblock+1,iblock+1));


        std::shared_ptr<BlockEpetra_Map> localRangeMap(
                                  new BlockEpetra_Map(localRangeBlockMaps));
        std::shared_ptr<BlockEpetra_Map> localDomainMap(
                                  new BlockEpetra_Map(localDomainBlockMaps));

        newPrec->setRangeMap(localRangeMap);
        newPrec->setDomainMap(localDomainMap);

        newPrec->updateApproximatedMomentumOperator();
        newPrec->updateApproximatedPressureMassOperator();
    }
    return newPrec;
}


void
GlobalSIMPLEOperator::
computeAm1BT(unsigned int rowIndex, unsigned int colIndex)
{
    ASSERT_PRE(colIndex >= M_nPrimalBlocks * 2, "Wrong col index!");
    ASSERT_PRE(rowIndex <  M_nPrimalBlocks * 2, "Wrong row index!");

    LifeChrono chrono;
    chrono.start();

    std::string msg = "Compute AM1BT for blocks (" + std::to_string(rowIndex) +
                      "," + std::to_string(colIndex) + ") ...";
    RedMA::printlog(RedMA::YELLOW, msg, M_verbose);

    matrixEpetraPtr_Type B = M_matrix.block(colIndex,rowIndex);
    matrixEpetraPtr_Type BT = M_matrix.block(rowIndex,colIndex);

    int numCols = BT->domainMap().mapSize();
    int numRows = BT->rangeMap().mapSize();

    // here we compute A^(-1) BT
    matrixEpetraPtr_Type resMatrix(new matrixEpetra_Type(BT->rangeMap(), numCols));
    resMatrix->zero();

    // auxiliary vector that we use to select the columns from BT
    vectorEpetra_Type aux(BT->domainMap());
    vectorEpetra_Type col(BT->rangeMap());
    vectorEpetra_Type res(BT->rangeMap());
    vectorEpetra_Type fakePressure(
                              M_matrix.block(rowIndex,rowIndex+1)->domainMap());
    fakePressure.zero();
    vectorEpetra_Type fakePressureResult(
                              M_matrix.block(rowIndex,rowIndex+1)->domainMap());
    fakePressureResult.zero();
    map_Type rangeMapRaw = col.epetraMap();
    unsigned int numElements = rangeMapRaw.NumMyElements();
    double dropTolerance = 1e-12;
    // retrieve vectors from matrix and apply simple operator to them
    for (unsigned int i = 0; i < numCols; i++)
    {
        aux.zero();
        col.zero();
        // do this only if the current processor owns the dof
        if (BT->domainMap().isOwned(i))
            aux[i] = 1.0;

        col = (*BT) * aux;

        // this sometimes fails for some reason. To debug.
        if (false)
        {
            mapEpetra_Type monolithicMap = BT->rangeMap();
            monolithicMap += M_matrix.block(rowIndex,rowIndex+1)->domainMap();

            vectorEpetra_Type subX(monolithicMap);
            vectorEpetra_Type subY(monolithicMap);

            subX.subset(col, monolithicMap, 0, 0);

            RedMA::CoutRedirecter ct;
            ct.redirect();

            M_invOperators[rowIndex/2]->ApplyInverse(subX.epetraVector(),
                                                     subY.epetraVector());

            ct.restore();

            res.subset(subY, monolithicMap, 0, 0);
        }
        else
            M_SingleOperators[rowIndex/2]->ApplyInverse(col, fakePressure,
                                                        res, fakePressureResult);

        for (unsigned int dof = 0; dof < numElements; dof++)
        {
            unsigned int gdof = rangeMapRaw.GID(dof);
            if (col.isGlobalIDPresent(gdof))
            {
                double value(col[gdof]);
                if (std::abs(value) > dropTolerance)
                {
                    resMatrix->addToCoefficient(gdof, i, value);
                }
            }
        }

    }
    M_comm->Barrier();
    std::shared_ptr<MapEpetra> domainMapPtr(new MapEpetra(BT->domainMap()));
    std::shared_ptr<MapEpetra> rangeMapPtr(new MapEpetra(BT->rangeMap()));

    resMatrix->globalAssemble(domainMapPtr, rangeMapPtr);

    M_Am1BT.block(rowIndex, colIndex) = resMatrix;

    msg = " done, in ";
    msg += std::to_string(chrono.diff());
    msg += " seconds\n";
    RedMA::printlog(RedMA::YELLOW, msg , M_verbose);
}

void
GlobalSIMPLEOperator::
computeGlobalSchurComplement()
{
    using namespace LifeV::MatrixEpetraStructuredUtility;
    typedef LifeV::MatrixEpetraStructuredView<double> MatrixView;

    std::string msg = "Assembling global Schur operator ...\n";
    LifeChrono chrono;
    chrono.start();
    RedMA::printlog(RedMA::GREEN, msg, M_verbose);

    M_Am1BT.resize(M_nBlockRows,M_nBlockCols);

    for (unsigned int i = 0; i < 2 * M_nPrimalBlocks; i++)
    {
        if (i % 2 == 0)
        {
            // loop on interfaces
            for (unsigned int j = M_nPrimalBlocks * 2; j < M_nBlockCols; j++)
            {
                if (M_oper(i,j) != nullptr && M_oper(i,i) != nullptr)
                {
                    computeAm1BT(i,j);
                }
            }
        }
    }


    M_globalSchurComplement.reset(new matrixEpetra_Type(*M_dualMap));

    LifeV::MatrixBlockStructure structure;
    structure.setBlockStructure(M_dimensionsInterfaces,
                                M_dimensionsInterfaces);


    unsigned int nDualBlocks = M_nBlockRows - 2 * M_nPrimalBlocks;

    for (unsigned int i = 0; i < nDualBlocks; i++)
    {
        for (unsigned int j = 0; j < nDualBlocks; j++)
        {
            matrixEpetraPtr_Type curMatrix(new matrixEpetra_Type(*M_dualMaps[i]));
            curMatrix->zero();
            for (unsigned int k = 0; k < 2 * M_nPrimalBlocks; k++)
            {
                matrixEpetraPtr_Type prod(new matrixEpetra_Type(*M_dualMaps[i]));
                prod->zero();
                if (M_matrix.block(i + M_nPrimalBlocks * 2, k) != nullptr &&
                    M_Am1BT.block(k, j + 2 * M_nPrimalBlocks)  != nullptr)
                {
                    M_matrix.block(i + M_nPrimalBlocks * 2, k)
                            ->multiply(false,
                                       *M_Am1BT.block(k, j + 2 * M_nPrimalBlocks),
                                       false,
                                       *prod,
                                       false);
                }
                prod->globalAssemble(M_dualMaps[j], M_dualMaps[i]);
                *curMatrix += *prod;
            }
            curMatrix->globalAssemble(M_dualMaps[j], M_dualMaps[i]);

            // copy into the global (small) matrix
            std::shared_ptr<MatrixView> blockGlobalView;
            blockGlobalView = createBlockView(M_globalSchurComplement,
                                              structure, i, j);

            LifeV::MatrixBlockStructure blockStructure;
            std::vector<unsigned int> rows(1), cols(1);
            rows[0] = M_dimensionsInterfaces[i];
            cols[0] = M_dimensionsInterfaces[j];
            blockStructure.setBlockStructure(rows, cols);

            std::shared_ptr<MatrixView> blockLocalView;
            blockLocalView = createBlockView(curMatrix, blockStructure, 0, 0);
            copyBlock(blockLocalView, blockGlobalView);
        }
    }
    M_globalSchurComplement->globalAssemble();
    // M_globalSchurComplement *= (-1);
    // create invertible approximated matrix
    M_approximatedGlobalSchurInverse.reset(new ApproximatedInvertibleMatrix);

    std::shared_ptr<Teuchos::ParameterList> globalSchurOptions;
    globalSchurOptions.reset(new
        Teuchos::ParameterList(M_solversOptions.sublist("GlobalSchurOperator")));
    M_approximatedGlobalSchurInverse->SetRowMatrix(
                                         M_globalSchurComplement->matrixPtr());
    M_approximatedGlobalSchurInverse->SetParameterList(*globalSchurOptions);
    M_approximatedGlobalSchurInverse->Compute();
    msg = "done, in ";
    msg += std::to_string(chrono.diff());
    msg += " seconds\n";
    RedMA::printlog(RedMA::GREEN, msg , M_verbose);
}

void
GlobalSIMPLEOperator::
fillComplete()
{
    // compute monolithic map
    M_monolithicMap.reset(new mapEpetra_Type());
    M_primalMap.reset(new mapEpetra_Type());
    M_dualMap.reset(new mapEpetra_Type());
    for (unsigned int iblock = 0; iblock < M_nBlockRows; iblock++)
    {
        for (unsigned int jblock = 0; jblock < M_nBlockCols; jblock++)
        {
            if (M_matrix.block(iblock,jblock) != nullptr)
            {
                M_allMaps.push_back(M_matrix.rangeMap(iblock,jblock));
                *M_monolithicMap +=  M_matrix.block(iblock,jblock)->rangeMap();
                if (iblock < 2 * M_nPrimalBlocks)
                {
                    M_primalMaps.push_back(M_matrix.rangeMap(iblock,jblock));
                    *M_primalMap += M_matrix.block(iblock,jblock)->rangeMap();
                }
                else
                {
                    M_dualMaps.push_back(M_matrix.rangeMap(iblock,jblock));
                    *M_dualMap += M_matrix.block(iblock,jblock)->rangeMap();
                    M_dimensionsInterfaces.push_back(
                        M_matrix.block(iblock,jblock)->rangeMap().mapSize());
                }
                break;
            }
        }
    }
}

void
GlobalSIMPLEOperator::
solveEveryPrimalBlock(const vectorEpetra_Type& X, vectorEpetra_Type &Y) const
{
    unsigned int offset = 0;
    for (unsigned int i = 0; i < M_nPrimalBlocks; i++)
    {
        mapEpetra_Type rangeVelocity = *M_matrix.rangeMap(i*2,i*2);
        mapEpetra_Type rangePressure;
        if (M_matrix.block(i*2+1,i*2) != nullptr)
            rangePressure = *M_matrix.rangeMap(i*2+1,i*2);
        else
            rangePressure = *M_matrix.rangeMap(i*2+1,i*2+1);

        if ((M_solvePrimalBlocksExactly && M_countIterations < M_numIterationExactSolve) ||
            (M_solvePrimalBlocksExactly && M_numIterationExactSolve < 0))
        {
            mapEpetra_Type monolithicMap = rangeVelocity;
            monolithicMap += rangePressure;

            vectorEpetra_Type subX(monolithicMap);
            vectorEpetra_Type subY(monolithicMap);

            subX.subset(X, monolithicMap, offset, 0);

            RedMA::CoutRedirecter ct;
            ct.redirect();

            M_invOperators[i]->ApplyInverse(subX.epetraVector(),
                                            subY.epetraVector());

            ct.restore();

            Y.subset(subY, monolithicMap, 0, offset);

            offset += monolithicMap.mapSize();
        }
        else
        {

        vectorEpetra_Type subVelocity(rangeVelocity);
        vectorEpetra_Type subPressure(rangePressure);
        subVelocity.zero();
        subPressure.zero();

        subVelocity.subset(X, rangeVelocity, offset, 0);
        subPressure.subset(X, rangePressure,
                           offset + rangeVelocity.mapSize(), 0);

        vectorEpetra_Type resVelocity(rangeVelocity);
        vectorEpetra_Type resPressure(rangePressure);

        M_SingleOperators[i]->ApplyInverse(subVelocity, subPressure,
                                           resVelocity, resPressure);

        Y.subset(resVelocity, rangeVelocity, 0, offset);
        Y.subset(resPressure, rangePressure, 0,
                 offset + rangeVelocity.mapSize());

        offset += rangeVelocity.mapSize() + rangePressure.mapSize();
        }
    }
}

void
GlobalSIMPLEOperator::
applyEveryB(const vectorEpetra_Type& X, vectorEpetra_Type &Y) const
{
    unsigned int offset = 0;

    std::vector<vectorEpetra_Type> Xs;
    for (unsigned int i = 0; i < M_nPrimalBlocks * 2; i++)
    {
        vectorEpetra_Type newVec(*M_allMaps[i]);
        newVec.subset(X, *M_allMaps[i], offset, 0);
        Xs.push_back(newVec);
        offset += M_allMaps[i]->mapSize();
    }

    offset = 0;
    for (unsigned int i = M_nPrimalBlocks * 2; i < M_nBlockRows; i++)
    {
        mapEpetra_Type curRange = *M_allMaps[i];
        vectorEpetra_Type subRes(curRange, LifeV::Unique);
        subRes.zero();
        for (unsigned int j = 0; j < M_nPrimalBlocks * 2; j++)
        {
            if (M_matrix.block(i,j))
            {
                subRes += (*M_matrix.block(i,j)) * Xs[j];
            }
        }
        Y.subset(subRes, curRange, 0, offset);
        offset += curRange.mapSize();
    }
}

void
GlobalSIMPLEOperator::
applyEveryBT(const vectorEpetra_Type& X, vectorEpetra_Type &Y) const
{
    unsigned int offset = 0;

    std::vector<vectorEpetra_Type> Xs;
    for (unsigned int i = M_nPrimalBlocks * 2; i < M_nBlockCols; i++)
    {
        vectorEpetra_Type newVec(*M_allMaps[i]);
        newVec.subset(X, *M_allMaps[i], offset, 0);
        Xs.push_back(newVec);
        offset += M_allMaps[i]->mapSize();
    }

    offset = 0;
    for (unsigned int i = 0; i < 2 * M_nPrimalBlocks; i++)
    {
        mapEpetra_Type curRange = *M_allMaps[i];
        vectorEpetra_Type subRes(curRange, LifeV::Unique);
        subRes.zero();
        for (unsigned int j = 2 * M_nPrimalBlocks; j < M_nBlockCols; j++)
        {
            if (M_matrix.block(i,j))
            {
                subRes += (*M_matrix.block(i,j)) * Xs[j-2*M_nPrimalBlocks];
            }
        }
        Y.subset(subRes, curRange, 0, offset);
        offset += curRange.mapSize();
    }
}

// application of the preconditioner corresponds to:
// primal part = Am1 x - A-1 BT(Sm1 B Am1 x - Sm1 y)
// dual part   = Sm1 (B Am1 x - y)
// Sm1 is the inverse of the Schur complement. We call z = Am1 x
int
GlobalSIMPLEOperator::
ApplyInverse(const vector_Type& X, vector_Type& Y) const
{
    ASSERT_PRE(X.NumVectors() == Y.NumVectors(),
               "X and Y must have the same number of vectors");
    if (M_nPrimalBlocks > 1)
    {
        const vectorEpetra_Type X_vectorEpetra(X, M_monolithicMap, LifeV::Unique);
        vectorEpetra_Type Y_vectorEpetra(M_monolithicMap, LifeV::Unique);

        vectorEpetra_Type X_primal(M_primalMap, LifeV::Unique);
        vectorEpetra_Type X_dual(M_dualMap, LifeV::Unique);
        X_primal.subset(X_vectorEpetra, *M_primalMap, 0, 0);
        X_dual.subset(X_vectorEpetra, *M_dualMap, M_primalMap->mapSize(), 0);

        // here we store the result
        vectorEpetra_Type Y_primal(M_primalMap, LifeV::Unique);
        vectorEpetra_Type Y_dual(M_dualMap, LifeV::Unique);

        vectorEpetra_Type Z(M_primalMap, LifeV::Unique);
        solveEveryPrimalBlock(X_primal, Z);

        vectorEpetra_Type Bz(M_dualMap, LifeV::Unique);
        applyEveryB(Z, Bz);

        M_approximatedGlobalSchurInverse->ApplyInverse((Bz - X_dual).epetraVector(),
                                                       Y_dual.epetraVector());

        vectorEpetra_Type BTy(M_primalMap, LifeV::Unique);
        applyEveryBT(Y_dual, BTy);

        vectorEpetra_Type Am1BTy(M_primalMap, LifeV::Unique);
        solveEveryPrimalBlock(BTy, Am1BTy);

        Y_primal = Z - Am1BTy;

        Y_vectorEpetra.subset(Y_primal, *M_primalMap, 0, 0);
        Y_vectorEpetra.subset(Y_dual, *M_dualMap, 0, M_primalMap->mapSize());
        Y = dynamic_cast<Epetra_MultiVector&>(Y_vectorEpetra.epetraVector());
    }
    else
    {
        M_SingleOperators[0]->ApplyInverse(X, Y);
    }
    M_countIterations++;
    return 0;
}
}
}
