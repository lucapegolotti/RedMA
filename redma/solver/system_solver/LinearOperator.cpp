#include "LinearOperator.hpp"

namespace RedMA
{

LinearOperator::
LinearOperator(const BM& matrix, SHP(BlockMaps) maps) :
  M_matrix(matrix),
  M_maps(maps)
{
    M_rangeMap.reset(new LifeV::BlockEpetra_Map(M_maps->M_rangeEpetraMaps));
    M_domainMap.reset(new LifeV::BlockEpetra_Map(M_maps->M_domainEpetraMaps));

    M_collapsedMatrix = M_maps->M_collapsedMatrix;
}

int
LinearOperator::
Apply(const super::vector_Type& X, super::vector_Type& Y) const
{
    using namespace LifeV;
    const std::unique_ptr<BlockEpetra_MultiVector> Xview(createBlockView(X, *M_domainMap));
    const std::unique_ptr<BlockEpetra_MultiVector> Yview(createBlockView(Y, *M_rangeMap));
    BlockEpetra_MultiVector tmpY(*M_rangeMap, X.NumVectors(), true);
    Yview->PutScalar(0.0);

    for (unsigned int i = 0; i < M_collapsedMatrix->nRows(); i++)
    {
        for (unsigned int j = 0; j < M_collapsedMatrix->nCols(); j++)
        {
            if (!M_collapsedMatrix->block(i,j)->isZero())
            {
                auto curBlock = std::static_pointer_cast<MATRIXEPETRA>(M_collapsedMatrix->block(i,j)->data());
                // EPETRA_CHK_ERR(curBlock->matrixPtr()->SetUseTranspose(true));
                // curBlock->spy("block" + std::to_string(i) + "_" + std::to_string(j));
                // std::cout << M_domainMap[j]->mapSize() << std::endl << std::flush;
                EPETRA_CHK_ERR(curBlock->matrixPtr()->Apply(Xview->block(j), tmpY.block(i)));
                EPETRA_CHK_ERR(Yview->block(i).Update(1.0, tmpY.block(i), 1.0));
                // EPETRA_CHK_ERR(curBlock->matrixPtr()->SetUseTranspose(false));
            }
        }
    }
    return 0;
}

}
