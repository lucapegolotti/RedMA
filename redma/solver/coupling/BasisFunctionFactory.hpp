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

#include <redma/geometry/BuildingBlock.hpp>

#include <redma/solver/coupling/BasisFunctionFunctor.hpp>
#include <redma/solver/coupling/DummyBasisFunction.hpp>
#include <redma/solver/coupling/ChebyshevBasisFunction.hpp>
#include <redma/solver/coupling/FourierBasisFunction.hpp>
#include <redma/solver/coupling/ZernikeBasisFunction.hpp>

#include <redma/utils/Exception.hpp>

namespace RedMA
{

SHP(BasisFunctionFunctor)
BasisFunctionFactory(const GetPot& datafile, GeometricFace inlet);

}  // namespace RedMA

#endif  // BASISFUNCTIONFACTORY_HPP
