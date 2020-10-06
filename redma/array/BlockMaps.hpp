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

#ifndef DIMENSION_HPP
#define DIMENSION_HPP

#include <redma/RedMA.hpp>
#include <redma/utils/Exception.hpp>
#include <redma/array/SparseMatrix.hpp>
#include <redma/array/BlockMatrix.hpp>

#include <lifev/core/array/MatrixEpetraStructured.hpp>
#include <lifev/core/array/VectorEpetraStructured.hpp>
#include <lifev/core/array/MatrixEpetraStructuredUtility.hpp>

namespace RedMA
{
//
// enum DimensionType{MAP, SIZE, MIXED, BLOCKDIMENSION};
//
// class Dimension
// {
// public:
//     Dimension(DimensionType type) : M_primaryType(type), M_secondaryType(type) {}
//
//     inline DimensionType primaryType() {return M_primaryType;}
//
//     inline DimensionType secondaryType() {return M_secondaryType;}
//
//     inline void setSecondaryType(DimensionType type) {M_secondaryType = type;}
//
//     inline unsigned int getNumElements() {return M_nElements;}
//
//     inline void setNumElements(unsigned int elements) {M_nElements = elements;}
//
//     inline SHP(MAPEPETRA) getMapEpetra() {return M_mapEpetra;}
//
//     inline void setMapEpetra(SHP(MAPEPETRA) mapEpetra) {M_mapEpetra = mapEpetra;}
//
//     inline SHP(Epetra_Map) getEpetraMap() {return M_epetraMap;}
//
//     inline SHP(Dimension)& dimension(unsigned int index) {throw new Exception("dimension() called on a simple Dimension object!");}
//
// protected:
//     // primary type = secondary type for single dimensions
//     // for blockDimensions, secondary type is the type of all the children (mixed if it varies)
//     DimensionType            M_primaryType;
//     DimensionType            M_secondaryType;
// private:
//     unsigned int             M_nElements;
//     SHP(MAPEPETRA)           M_mapEpetra;
//     SHP(Epetra_Map)          M_epetraMap;
// };

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

// SHP(aMatrix) collapseBlocks(const SHP(aMatrix)& matrix,
//                             const BlockDimension& maps);

SHP(BlockMatrix) collapseBlocks(SHP(BlockMatrix) matrix,
                                std::vector<unsigned int>& dimensionsRowsBlock,
                                std::vector<unsigned int>& dimensionsColsBlock);

SHP(SparseMatrix) blockMatrixToSparseMatrix(SHP(BlockMatrix) matrix);

SHP(SparseMatrix) denseMatrixToSparseMatrix(SHP(DenseMatrix) matrix);

// void buildBlockMaps(std::vector<SHP(MAPEPETRA)> rangeMaps,
//                     std::vector<SHP(MAPEPETRA)> domainMaps,
//                     SHP(MAPEPETRA) monolithicRangeMaps,
//                     SHP(MAPEPETRA) monolithicDomainMaps);

SHP(VECTOREPETRA) getEpetraVector(const SHP(aVector)& vector,
                                  const BlockMaps& maps);
//
SHP(aVector) getBlockVector(const SHP(VECTOREPETRA)& vector,
                            const BlockMaps& maps);

}

#endif // BLOCKSIMENSIONS_HPP
