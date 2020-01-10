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
#include <AbstractMatrix.hpp>
#include <Vector.hpp>

#include <lifev/core/array/MatrixEpetra.hpp>

namespace RedMA
{

template<class InMatrixType>
class Matrix : public AbstractMatrix
{
    typedef std::shared_ptr<InMatrixType>            InMatrixTypePtr;
public:
    Matrix();

    // soft copy
    Matrix(const InMatrixTypePtr& matrix);

    // soft copy
    Matrix(const Matrix<InMatrixType>& other);

    // hard copy
    Matrix<InMatrixType>& operator=(const Matrix<InMatrixType>& other);

    // hard copy
    Matrix<InMatrixType>& operator=(const InMatrixTypePtr& matrix);

    Matrix<InMatrixType>& operator+=(const Matrix<InMatrixType>& other);

    Matrix<InMatrixType>& operator*=(const double& coeff);

    template<class InVectorType>
    Vector<InVectorType> operator*(const Vector<InVectorType>& vector);

    InMatrixTypePtr& get();

    InMatrixTypePtr get() const;

    bool isZero() const;

    bool isNull() const;

private:
    InMatrixTypePtr                     M_inMatrix;
};

}  // namespace RedMA

#include <Matrix_imp.hpp>

#endif  // MATRIXEP_HPP
