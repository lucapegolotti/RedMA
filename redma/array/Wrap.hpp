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

#ifndef WRAP_HPP
#define WRAP_HPP

#include <redma/array/DenseMatrix.hpp>
#include <redma/array/DenseVector.hpp>
#include <redma/array/SparseMatrix.hpp>
#include <redma/array/DistributedVector.hpp>

namespace RedMA
{

// Wrap a raw dense matrix into its wrapper.
shp<DenseMatrix> wrap(shp<DENSEMATRIX> matrix);

// Wrap a raw dense vector into its wrapper.
shp<DenseVector> wrap(shp<DENSEVECTOR> vector);

// Wrap a raw sparse matrix into its wrapper.
shp<SparseMatrix> wrap(shp<MATRIXEPETRA> matrix);

// Wrap a raw distributed matrix into its wrapper.
shp<DistributedVector> wrap(shp<VECTOREPETRA> vector);

}

#endif // WRAP_HPP
