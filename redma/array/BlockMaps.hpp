// Reduced Modeling of Arteries (RedMA)
// Copyright (C) 2019  Luca Pegolotti
//
// RedMA is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// RedMA is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with this program.  If not, see <http://www.gnu.org/licenses/>.

#ifndef BLOCKMAPS_HPP
#define BLOCKMAPS_HPP

#include <redma/RedMA.hpp>
#include <redma/utils/Exception.hpp>
#include <redma/array/SparseMatrix.hpp>
#include <redma/array/BlockMatrix.hpp>

#include <lifev/core/array/MatrixEpetraStructured.hpp>
#include <lifev/core/array/VectorEpetraStructured.hpp>
#include <lifev/core/array/MatrixEpetraStructuredUtility.hpp>

namespace RedMA
{

class BlockMaps
{
public:
    BlockMaps() {}

    BlockMaps(shp<BlockMatrix> matrix) {createFromBlockMatrix(matrix);}

    void createFromBlockMatrix(shp<BlockMatrix> matrix);

    void updateCollapsedMatrix(shp<BlockMatrix> matrix);

    shp<BlockMatrix>                 M_collapsedMatrix;
    std::vector<shp<MAPEPETRA>>      M_rangeMaps;
    std::vector<shp<MAPEPETRA>>      M_domainMaps;
    std::vector<shp<Epetra_Map>>     M_rangeEpetraMaps;
    std::vector<shp<Epetra_Map>>     M_domainEpetraMaps;
    shp<MAPEPETRA>                   M_monolithicRangeMap;
    shp<MAPEPETRA>                   M_monolithicDomainMap;
    std::vector<unsigned int>        M_dimensionsRows;
    std::vector<unsigned int>        M_dimensionsCols;
    std::vector<unsigned int>        M_dimensionsRowsBlock;
    std::vector<unsigned int>        M_dimensionsColsBlock;
};

shp<BlockMatrix> collapseBlocks(shp<BlockMatrix> matrix,
                                std::vector<unsigned int>& dimensionsRowsBlock,
                                std::vector<unsigned int>& dimensionsColsBlock);

shp<SparseMatrix> blockMatrixToSparseMatrix(shp<BlockMatrix> matrix);

shp<SparseMatrix> denseMatrixToSparseMatrix(shp<DenseMatrix> matrix);

shp<VECTOREPETRA> getEpetraVector(const shp<aVector>& vector,
                                  const BlockMaps& maps);

shp<aVector> getBlockVector(const shp<VECTOREPETRA>& vector,
                            const BlockMaps& maps);
}

#endif // BLOCKMAPS_HPP
