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

#ifndef STOKESMODEL_HPP
#define STOKESMODEL_HPP

#include <redma/RedMA.hpp>
#include <redma/assemblers/abstract/aAssembler.hpp>
#include <redma/array/DistributedVector.hpp>
#include <redma/array/SparseMatrix.hpp>
#include <redma/geometry/BuildingBlock.hpp>

// #include <redma/reduced_basis/BlockMDEIM.hpp>

#include <lifev/eta/expression/Integrate.hpp>
#include <lifev/core/filter/Exporter.hpp>
#include <lifev/core/filter/ExporterVTK.hpp>
#include <lifev/core/filter/ExporterHDF5.hpp>

namespace RedMA
{

class StokesModel
{
public:
    StokesModel(const DataContainer& data, shp<TreeNode> treeNode);

protected:

    shp<BlockVector> buildZeroVector() const;

};

}

#endif // STOKESMODEL_HPP