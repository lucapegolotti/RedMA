#include "BlockMaps.hpp"

namespace RedMA
{

template <>
void
BlockMaps<MatrixEp>::
generateMaps()
{
    unsigned int nrows = M_matrix.nRows();
    unsigned int ncols = M_matrix.nCols();

    M_rangeMaps.resize(nrows);
    M_domainMaps.resize(ncols);

    for (unsigned int i = 0; i < nrows; i++)
    {
        for (unsigned int j = 0; j < ncols; j++)
        {
            if (!M_matrix.block(i,j).isNull())
            {
                M_rangeMaps[i] = M_matrix.block(i,j).data()->rangeMap().map(LifeV::Unique);
                M_domainMaps[i] = M_matrix.block(i,j).data()->domainMap().map(LifeV::Unique);
            }
        }
    }
}

template <>
void
BlockMaps<BlockMatrix<MatrixEp>>::
generateMaps()
{
    unsigned int nrows = M_matrix.nRows();
    unsigned int ncols = M_matrix.nCols();

    // first we deal with ranges
    for (unsigned int i = 0; i < nrows; i++)
    {
        for (unsigned int j = 0; j < ncols; j++)
        {
            if (!M_matrix.block(i,j).isNull())
            {
                BlockMaps<MatrixEp> localMaps(M_matrix.block(i,j));
                auto ranges = localMaps.getRangeMaps();

                for (auto map : ranges)
                    M_rangeMaps.push_back(map);
                break;
            }
        }
    }

    for (unsigned int j = 0; j < ncols; j++)
    {
        for (unsigned int i = 0; i < nrows; i++)
        {
            if (!M_matrix.block(i,j).isNull())
            {
                BlockMaps<MatrixEp > localMaps(M_matrix.block(i,j));
                auto ranges = localMaps.getDomainMaps();

                for (auto map : ranges)
                    M_domainMaps.push_back(map);
                break;
            }
        }
    }
}

BlockMatrix<MatrixEp> collapseBlocks(const BlockMatrix<BlockMatrix<MatrixEp>>& matrix,
                                     const BlockMaps<BlockMatrix<MatrixEp>>& maps)
{
    if (!matrix.isFinalized())
    {
        throw new Exception("Matrix must be finalized in order to collapse blocks!");
    }

    auto rangeMaps = maps.getRangeMaps();
    auto domainMaps = maps.getRangeMaps();

    unsigned int nrows = rangeMaps.size();
    unsigned int ncols = domainMaps.size();

    BlockMatrix<MatrixEp> retMatrix;
    retMatrix.resize(nrows, ncols);

    unsigned int offsetrow = 0;
    for (unsigned int i = 0; i < matrix.nRows(); i++)
    {
        unsigned int offsetcol = 0;
        unsigned int curnrows = matrix.block(i,0).nRows();
        for (unsigned int j = 0; j < matrix.nCols(); j++)
        {
            auto cursubMatrix = matrix.block(i,j);

            for (unsigned int ii = 0; ii < curnrows; ii++)
            {
                for (unsigned int jj = 0; jj < cursubMatrix.nCols(); jj++)
                {
                    retMatrix.block(offsetrow+ii,offsetcol+jj).softCopy(cursubMatrix.block(ii,jj));
                }
            }
            offsetcol += cursubMatrix.nCols();
        }
        offsetrow += curnrows;
    }

    return retMatrix;
}

}
