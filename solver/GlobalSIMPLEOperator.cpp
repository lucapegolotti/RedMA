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
GlobalSIMPLEOperator() :
  M_label("GlobalSIMPLEOperator"),
  M_useTranspose(false)
{

}

GlobalSIMPLEOperator::~GlobalSIMPLEOperator()
{

}

void
GlobalSIMPLEOperator::
setUp(RedMA::GlobalBlockMatrix matrix,
      const commPtr_Type & comm)
{
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

            localRangeBlockMaps[1].reset(new
                    Epetra_Map(blockOper(iblock,iblock+1)->OperatorRangeMap()));

            localDomainBlockMaps[0].reset(new
                    Epetra_Map(blockOper(iblock,iblock)->OperatorDomainMap()));

            localDomainBlockMaps[1].reset(new
                    Epetra_Map(blockOper(iblock+1,iblock)->OperatorDomainMap()));

            // here we rely on the structure of the global matrix
            PreconditionerPtr newPrec;
            newPrec.reset(Operators::NSPreconditionerFactory::
                          instance().createObject("SIMPLE"));
            newPrec->setOptions(M_solversOptions);
            newPrec->setUp(matrix.block(iblock,iblock),
                           matrix.block(iblock+1,iblock),
                           matrix.block(iblock,iblock+1));

            std::shared_ptr<BlockEpetra_Map> localRangeMap(
                                      new BlockEpetra_Map(localRangeBlockMaps));
            std::shared_ptr<BlockEpetra_Map> localDomainMap(
                                      new BlockEpetra_Map(localDomainBlockMaps));

            newPrec->setRangeMap(localRangeMap);
            newPrec->setDomainMap(localDomainMap);
            newPrec->updateApproximatedMomentumOperator();
            newPrec->updateApproximatedSchurComplementOperator();
            M_SIMPLEOperators.push_back(newPrec);

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

    M_approximatedBAm1Binverses.resize(M_nBlockRows, M_nBlockCols);
    computeBAm1BT_inverse(0,4);
}

void
GlobalSIMPLEOperator::
computeBAm1BT_inverse(unsigned int rowIndex, unsigned int colIndex)
{
    ASSERT_PRE(colIndex >= M_nPrimalBlocks * 2, "Wrong col index!");
    ASSERT_PRE(rowIndex <  M_nPrimalBlocks * 2, "Wrong row index!");

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

        M_SIMPLEOperators[rowIndex / 2]->ApplyInverse(col, fakePressure,
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

    // compute product B * resMatrix
    matrixEpetraPtr_Type BtildeBT(new matrixEpetra_Type(BT->domainMap()));
    B->multiply(false,* resMatrix, false, *BtildeBT, false);

    BtildeBT->globalAssemble(std::make_shared<MapEpetra>(BT->domainMap()),
                             std::make_shared<MapEpetra>(BT->domainMap()));

    // create invertible approximated matrix
    M_approximatedBAm1Binverses(rowIndex, colIndex).reset(
                                new Operators::ApproximatedInvertibleRowMatrix);

    std::shared_ptr<Teuchos::ParameterList> globalSchurOptions;
    globalSchurOptions.reset(new
        Teuchos::ParameterList(M_solversOptions.sublist("GlobalSchurOperator")));
    M_approximatedBAm1Binverses(rowIndex, colIndex)->
                                            SetRowMatrix(BtildeBT->matrixPtr());
    M_approximatedBAm1Binverses(rowIndex, colIndex)->
                                          SetParameterList(*globalSchurOptions);
    M_approximatedBAm1Binverses(rowIndex, colIndex)->Compute();
}

void
GlobalSIMPLEOperator::
fillComplete()
{
    // Filling the empty blocks with null operators
    for (UInt iblock = 0; iblock < M_nBlockRows; ++iblock)
        for (UInt jblock = 0; jblock < M_nBlockCols; ++jblock)
        {
            if (M_oper(iblock,jblock).get() == 0)
            {
                NullOperator * nullOp(new NullOperator);
                nullOp->setUp(M_domainMap->blockMap(jblock), M_rangeMap->blockMap(iblock));
                M_oper(iblock,jblock).reset(nullOp);
            }
            ASSERT(M_rangeMap->blockMap(iblock)->
                   PointSameAs(M_oper(iblock,jblock)->OperatorRangeMap()),
                   "Wrong range map");
            ASSERT(M_domainMap->blockMap(jblock)->
                   PointSameAs(M_oper(iblock,jblock)->OperatorDomainMap()),
                   "Wrong domain map");
        }
}


int GlobalSIMPLEOperator::ApplyInverse(const vector_Type& X, vector_Type& Y) const
{
    ASSERT_PRE(X.NumVectors() == Y.NumVectors(), "X and Y must have the same number of vectors");

    Y = X;

    return 0;
}
}
}
