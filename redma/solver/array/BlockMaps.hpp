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
#include <redma/solver/array/BlockMatrix.hpp>
#include <redma/solver/array/MatrixEp.hpp>

namespace RedMA
{

template <class InMatrixType>
class BlockMaps
{
public:
    BlockMaps(const BlockMatrix<InMatrixType>& matrix);

    inline std::vector<SHP(Epetra_Map)> getRangeMaps() const {return M_rangeMaps;}

    inline std::vector<SHP(Epetra_Map)> getDomainMaps() const {return M_domainMaps;}

private:
    void generateMaps();

    BlockMatrix<InMatrixType>  M_matrix;
    std::vector<SHP(Epetra_Map)> M_rangeMaps;
    std::vector<SHP(Epetra_Map)> M_domainMaps;
};

BlockMatrix<MatrixEp> collapseBlocks(const BlockMatrix<BlockMatrix<MatrixEp>>& matrix,
                                     const BlockMaps<BlockMatrix<MatrixEp>>& maps);

}

#include "BlockMaps_imp.hpp"

#endif // BLOCKMAPS_HPP
