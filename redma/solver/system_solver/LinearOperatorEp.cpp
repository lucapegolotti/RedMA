#include "LinearOperatorEp.hpp"

namespace RedMA
{

LinearOperatorEp::
LinearOperatorEp()
{
}

void
LinearOperatorEp::
setup(BM matrix, EPETRACOMM comm)
{
    using namespace LifeV;

    M_matrix.softCopy(matrix);
    M_comm = comm;

    std::vector<super::mapPtr_Type> rangeBlockMaps(M_matrix.nRows());
    std::vector<super::mapPtr_Type> domainBlockMaps(M_matrix.nCols());

    for (unsigned int i = 0; i < M_matrix.nRows(); i++)
    {
        SHP(LifeV::MapEpetra) curMap;
        M_matrix.getRowProperty(curMap, i);

        rangeBlockMaps[i] = curMap->map(Unique);
    }

    for (unsigned int j = 0; j < M_matrix.nCols(); j++)
    {
        SHP(LifeV::MapEpetra) curMap;
        M_matrix.getColProperty(curMap, j);

        domainBlockMaps[j] = curMap->map(Unique);
    }

    M_rangeMap.reset(new BlockEpetra_Map(rangeBlockMaps));
    M_domainMap.reset(new BlockEpetra_Map(domainBlockMaps));
}

int
LinearOperatorEp::
Apply(const super::vector_Type& X, super::vector_Type& Y) const
{
    using namespace LifeV;
    const std::unique_ptr<BlockEpetra_MultiVector> Xview(createBlockView(X, *M_domainMap));
    const std::unique_ptr<BlockEpetra_MultiVector> Yview(createBlockView(Y, *M_rangeMap));
    BlockEpetra_MultiVector tmpY(*M_rangeMap, X.NumVectors(), true);

    for (unsigned int i = 0; i < M_matrix.nRows(); i++)
    {
        for (unsigned int j = 0; j < M_matrix.nCols(); j++)
        {
            auto& curmat = M_matrix.block(i,j).data()->matrixPtr();
            curmat->SetUseTranspose(true);
            curmat->Apply(Xview->block(j), tmpY.block(i));
            Yview->block(i).Update(1.0, tmpY.block(i), 1.0);
            curmat->SetUseTranspose(false);
        }
    }
}

}
