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

#ifndef MATRIX_HPP
#define MATRIX_HPP

#include <Exception.hpp>

#include <lifev/core/array/MapVector.hpp>
#include <lifev/core/array/VectorEpetraStructured.hpp>

#include <AbstractVector.hpp>

#include <boost/numeric/ublas/matrix.hpp>

namespace RedMA
{

// // this is simply the interface of a wrapper to wathever type of matrix we use
// template<class InMatrixType>
// class Matrix
// {
//     typedef std::shared_ptr<InMatrixType>                       InMatrixPtr;
//     typedef LifeV::VectorEpetra                                 VectorEpetra;
//     typedef std::shared_ptr<VectorEpetra>                       VectorEpetraPtr;
//     typedef LifeV::MapEpetra                                    Map;
//     typedef std::shared_ptr<Map>                                MapPtr;
//
// public:
//     Matrix();
//
//     virtual Matrix(const Matrix& other) = 0;
//
//     virtual void add(const Matrix& other) = 0;
//
//     virtual void VectorEpetra operator*(const Vector& vector) = 0;
//
//     Matrix& operator*=(const double& coeff) = 0;
//
//     void zero();
//
//     InMatrixPtr get();
//
// protected:
//     InMatrixPtr          M_matrix;
// };
//
// #include "Matrix_imp.hpp";

// this is simply the interface of a wrapper to wathever type of matrix we use
class AbstractMatrix
{

};

}  // namespace RedMA

#endif  // MATRIX_HPP
