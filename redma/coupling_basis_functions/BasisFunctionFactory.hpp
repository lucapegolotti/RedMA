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

#ifndef BASISFUNCTIONFACTORY_HPP
#define BASISFUNCTIONFACTORY_HPP

#include <memory>

#include <redma/RedMA.hpp>

#include <redma/geometry/building_blocks/BuildingBlock.hpp>

#include <redma/coupling_basis_functions/BasisFunctionFunctor.hpp>
#include <redma/coupling_basis_functions/DummyBasisFunction.hpp>
#include <redma/coupling_basis_functions/ChebyshevBasisFunction.hpp>
#include <redma/coupling_basis_functions/FourierBasisFunction.hpp>
#include <redma/coupling_basis_functions/ZernikeBasisFunction.hpp>
#include <redma/coupling_basis_functions/FourierRingBasisFunction.hpp>

#include <redma/utils/Exception.hpp>

namespace RedMA
{

/// \brief Factory for BasisFunctionFunctor.
shp<BasisFunctionFunctor>
BasisFunctionFactory(const GetPot& datafile, GeometricFace inlet,
                     bool isBoundary = false, bool isRing = false,
                     const double mesh_size = 0.0);

}  // namespace RedMA

#endif  // BASISFUNCTIONFACTORY_HPP
