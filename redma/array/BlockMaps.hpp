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

/// Class for handling the dimensions of blocks in block structures.
class BlockMaps
{
public:
    /// Default constructor.
    BlockMaps() {}

    /*! \brief Constructor.
     *
     * Internally, we call the method createFromBlockMatrix.
     * \param The input matrix.
     */
    BlockMaps(shp<BlockMatrix> matrix) {createFromBlockMatrix(matrix);}

    /*! \brief Constructs block maps from an input matrix.
     *
     * \param The input matrix.
     */
    void createFromBlockMatrix(shp<BlockMatrix> matrix);


    /*! \brief Updates the internal collapsed matrix.
     *
     * \param The input matrix.
     */
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

/* Given a block matrix, this method returns a "flattened" block matrix in which
 * each block is a non block structure. It also returns the vectors of dimensions
 * of rows and columns.
 */
shp<BlockMatrix> collapseBlocks(shp<BlockMatrix> matrix,
                                std::vector<unsigned int>& dimensionsRowsBlock,
                                std::vector<unsigned int>& dimensionsColsBlock);

// Converts a block matrix to a sparse one.
shp<SparseMatrix> blockMatrixToSparseMatrix(shp<BlockMatrix> matrix);

// Converts a dense matrix to a sparse one.
shp<SparseMatrix> denseMatrixToSparseMatrix(shp<DenseMatrix> matrix);

// Converts a block matrix to a dense one.
shp<DenseMatrix> blockMatrixToDenseMatrix(shp<BlockMatrix> matrix);

// Converts a block vector to a dense one.
shp<DenseVector> blockVectorToDenseVector(shp<BlockVector> vector);

/* Given an input block vector and the corresponding maps, this function returns
 * a monolithic epetra vector.
 */
shp<VECTOREPETRA> getEpetraVector(const shp<aVector>& vector,
                                  const BlockMaps& maps);

/* Given an input epetra vector and block maps, this function returns
 * a block vector with row dimensions corresponding to the input maps.
 */
shp<aVector> getBlockVector(const shp<VECTOREPETRA>& vector,
                            const BlockMaps& maps);
}

#endif // BLOCKMAPS_HPP
