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
setUp(operatorPtrContainer_Type blockOper,
      const commPtr_Type & comm)
{
    M_nBlockRows = blockOper.size1();
    M_nBlockCols = blockOper.size2();
    M_comm = comm;

    BlockEpetra_Map::mapPtrContainer_Type rangeBlockMaps(M_nBlockRows);
    BlockEpetra_Map::mapPtrContainer_Type domainBlockMaps(M_nBlockCols);

    for (UInt iblock=0; iblock < M_nBlockRows; ++iblock)
        for (UInt jblock=0; jblock < M_nBlockCols; ++jblock)
        {
            if (blockOper(iblock,jblock) != 0 && rangeBlockMaps[iblock]==0)
            {
                rangeBlockMaps[iblock].reset(new
                      Epetra_Map(blockOper(iblock,jblock)->OperatorRangeMap()));
                jblock = M_nBlockCols;
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
