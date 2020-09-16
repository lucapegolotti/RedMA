#include "LinearOperator.hpp"

namespace RedMA
{

LinearOperator::
LinearOperator(const BM& matrix) :
  M_matrix(matrix)
  // M_maps(new BlockMaps<BlockMatrix<MatrixEp>>(matrix))
{
    // M_rangeMap.reset(new LifeV::BlockEpetra_Map(M_maps->getRangeMaps()));
    // M_domainMap.reset(new LifeV::BlockEpetra_Map(M_maps->getDomainMaps()));

    // M_collapsedMatrix = collapseBlocks(M_matrix, *M_maps);
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
            // EPETRA_CHK_ERR(M_collapsedMatrix.block(i,j).data()->matrixPtr()->SetUseTranspose(true));
            // EPETRA_CHK_ERR(M_collapsedMatrix.block(i,j).data()->matrixPtr()->Apply(Xview->block(j), tmpY.block(i)));
            // EPETRA_CHK_ERR(Yview->block(i).Update(1.0, tmpY.block(i), 1.0));
            // EPETRA_CHK_ERR(M_collapsedMatrix.block(i,j).data()->matrixPtr()->SetUseTranspose(false));
        }
    }

    return 0;
}

}
