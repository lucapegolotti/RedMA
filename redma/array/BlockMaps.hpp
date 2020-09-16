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
#include <redma/array/BlockMatrix.hpp>
#include <redma/array/SparseMatrix.hpp>

#include <lifev/core/array/MatrixEpetraStructured.hpp>
#include <lifev/core/array/VectorEpetraStructured.hpp>
#include <lifev/core/array/MatrixEpetraStructuredUtility.hpp>

namespace RedMA
{

class BlockMaps
{
public:
    BlockMaps(const BlockMatrix& matrix) {};

    inline std::vector<SHP(Epetra_Map)> getRangeMaps() const {return M_rangeMaps;}

    inline std::vector<SHP(Epetra_Map)> getDomainMaps() const {return M_domainMaps;}

    inline std::vector<SHP(MAPEPETRA)> getRangeMapsEpetra() const {return M_rangeMapsEpetra;}

    inline std::vector<SHP(MAPEPETRA)> getDomainMapsEpetra() const {return M_domainMapsEpetra;}

    inline std::vector<unsigned int> getInnerRows() const {return M_nInnerRows;}

    inline std::vector<unsigned int> getInnerCols() const {return M_nInnerCols;}

    SHP(MAPEPETRA) getMonolithicRangeMapEpetra() const {};

    SHP(MAPEPETRA) getMonolithicDomainMapEpetra() const {};

private:
    void generateMaps() {};

    SHP(BlockMatrix)             M_matrix;
    std::vector<SHP(Epetra_Map)> M_rangeMaps;
    std::vector<SHP(Epetra_Map)> M_domainMaps;
    std::vector<SHP(MAPEPETRA)>  M_rangeMapsEpetra;
    std::vector<SHP(MAPEPETRA)>  M_domainMapsEpetra;

    std::vector<unsigned int>    M_nInnerRows;
    std::vector<unsigned int>    M_nInnerCols;
};

SHP(aMatrix) collapseBlocks(const SHP(aMatrix)& matrix,
                            const BlockMaps& maps);

SHP(VECTOREPETRA) getEpetraVector(const SHP(aVector)& vector,
                                  const BlockMaps& maps);

SHP(aVector) getBlockVector(const SHP(VECTOREPETRA)& vector,
                                const BlockMaps& maps);

}

#endif // BLOCKMAPS_HPP
