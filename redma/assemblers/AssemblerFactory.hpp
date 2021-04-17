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

#ifndef ASSEMBLERFACTORY_HPP
#define ASSEMBLERFACTORY_HPP

#include <redma/RedMA.hpp>
#include <redma/assemblers/abstract/aAssembler.hpp>
#include <redma/assemblers/finite_element/StokesAssemblerFE.hpp>
#include <redma/assemblers/reduced_basis/StokesAssemblerRB.hpp>
#include <redma/assemblers/finite_element/NavierStokesAssemblerFE.hpp>
#include <redma/assemblers/reduced_basis/NavierStokesAssemblerRB.hpp>

#include <redma/assemblers/finite_element/FSIAssemblerFE.hpp>
#include <redma/utils/Exception.hpp>
#include <redma/geometry/building_blocks/BuildingBlock.hpp>
#include <redma/problem/DataContainer.hpp>

namespace RedMA
{

/*! \brief Factory to create assemblers.
 *
 * \param data The DataContainer of the problem.
 * \param treeNode Shared pointer to the tree node.
 * \return Shared pointer to the abstract assembler.
 */
shp<aAssembler>
AssemblerFactory(const DataContainer& data,
                 shp<TreeNode> treeNode);

}

#endif // ASSEMBLERFACTORY_HPP
