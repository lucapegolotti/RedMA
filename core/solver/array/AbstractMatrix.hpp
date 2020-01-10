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

#ifndef ABSTRACTMATRIX_HPP
#define ABSTRACTMATRIX_HPP

#include <Exception.hpp>

#include <memory>

#include <lifev/core/array/MapVector.hpp>
#include <lifev/core/array/VectorEpetraStructured.hpp>

#include <AbstractVector.hpp>

#include <boost/numeric/ublas/matrix.hpp>

namespace RedMA
{

// this should be the interface of a wrapper to wathever type of matrix we use
class AbstractMatrix
{

};

}  // namespace RedMA

#endif  // ABSTRACTMATRIX_HPP
