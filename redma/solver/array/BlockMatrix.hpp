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

#ifndef BLOCKMATRIX_HPP
#define BLOCKMATRIX_HPP

#include <redma/utils/Exception.hpp>
#include <redma/solver/array/aMatrix.hpp>

#include <boost/numeric/ublas/matrix.hpp>

namespace RedMA
{

template <class InMatrixType>
class BlockMatrix : public aMatrix
{
    typedef boost::numeric::ublas::matrix<InMatrixType>       Grid;
public:
    BlockMatrix(const unsigned int& nRows, const unsigned int& nCols);

    virtual BlockMatrix<InMatrixType> operator+(const BlockMatrix<InMatrixType>& other) const;

    virtual BlockMatrix<InMatrixType>& operator+=(const BlockMatrix<InMatrixType>& other);

    virtual BlockMatrix<InMatrixType>& operator*=(const double& coeff);

    template <class VectorType>
    VectorType operator*(const VectorType& vector);

    void resize(const unsigned int& nRows, const unsigned int& nCols);

    void hardCopy(const BlockMatrix<InMatrixType>& other);

    InMatrixType& block(const unsigned int& iblock, const unsigned int& jblock);

    InMatrixType block(const unsigned int& iblock, const unsigned int& jblock) const;

private:
    Grid            M_matrixGrid;
    unsigned int    M_nRows;
    unsigned int    M_nCols;
};

}

#include "BlockMatrix_imp.hpp"

#endif // BLOCKMATRIX_HPP
