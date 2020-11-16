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
    BlockMaps();

    BlockMaps(SHP(BlockMatrix) matrix) {createFromBlockMatrix(matrix);}

    void createFromBlockMatrix(SHP(BlockMatrix) matrix);

    void updateCollapsedMatrix(SHP(BlockMatrix) matrix);

    SHP(BlockMatrix)                 M_collapsedMatrix;
    std::vector<SHP(MAPEPETRA)>      M_rangeMaps;
    std::vector<SHP(MAPEPETRA)>      M_domainMaps;
    std::vector<SHP(Epetra_Map)>     M_rangeEpetraMaps;
    std::vector<SHP(Epetra_Map)>     M_domainEpetraMaps;
    SHP(MAPEPETRA)                   M_monolithicRangeMap;
    SHP(MAPEPETRA)                   M_monolithicDomainMap;
    std::vector<unsigned int>        M_dimensionsRows;
    std::vector<unsigned int>        M_dimensionsCols;
    std::vector<unsigned int>        M_dimensionsRowsBlock;
    std::vector<unsigned int>        M_dimensionsColsBlock;
};

SHP(BlockMatrix) collapseBlocks(SHP(BlockMatrix) matrix,
                                std::vector<unsigned int>& dimensionsRowsBlock,
                                std::vector<unsigned int>& dimensionsColsBlock);

SHP(SparseMatrix) blockMatrixToSparseMatrix(SHP(BlockMatrix) matrix);

SHP(SparseMatrix) denseMatrixToSparseMatrix(SHP(DenseMatrix) matrix);

SHP(VECTOREPETRA) getEpetraVector(const SHP(aVector)& vector,
                                  const BlockMaps& maps);

SHP(aVector) getBlockVector(const SHP(VECTOREPETRA)& vector,
                            const BlockMaps& maps);
}

#endif // BLOCKMAPS_HPP
