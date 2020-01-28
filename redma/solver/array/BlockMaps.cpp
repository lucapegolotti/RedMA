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

}
