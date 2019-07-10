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
    operatorPtrContainer_Type blockOper = matrix.getGrid();

    M_nBlockRows = blockOper.size1();
    M_nBlockCols = blockOper.size2();
    M_comm = comm;

    BlockEpetra_Map::mapPtrContainer_Type rangeBlockMaps(M_nBlockRows);
    BlockEpetra_Map::mapPtrContainer_Type domainBlockMaps(M_nBlockCols);

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
                           matrix.block(iblock,iblock+1),
                           matrix.block(iblock+1,iblock));

            std::shared_ptr<BlockEpetra_Map> localRangeMap(
                                      new BlockEpetra_Map(localRangeBlockMaps));
            std::shared_ptr<BlockEpetra_Map> localDomainMap(
                                      new BlockEpetra_Map(localDomainBlockMaps));

            newPrec->setRangeMap(localRangeMap);
            newPrec->setDomainMap(localDomainMap);
            // newPrec->updateApproximatedMomentumOperator();
            // newPrec->updateApproximatedSchurComplementOperator();
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
