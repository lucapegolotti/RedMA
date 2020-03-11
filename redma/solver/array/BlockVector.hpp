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

#ifndef BLOCKVECTOR_HPP
#define BLOCKVECTOR_HPP

#include <redma/utils/Exception.hpp>
#include <redma/solver/array/DenseVector.hpp>

#include <boost/numeric/ublas/matrix.hpp>

#include <math.h>

namespace RedMA
{

template <class InVectorType>
class BlockVector
{
    typedef boost::numeric::ublas::matrix<InVectorType>       Grid;

public:
    typedef InVectorType                                      InnerType;

    BlockVector();

    BlockVector(const BlockVector& other);

    BlockVector(BlockVector& other);

    BlockVector(const unsigned int& nRows);

    BlockVector operator*(const double& coeff) const;

    BlockVector& operator*=(const double& coeff);

    BlockVector& operator+=(const BlockVector<InVectorType>& other);

    BlockVector& operator-=(const BlockVector<InVectorType>& other);

    BlockVector operator+(const BlockVector<InVectorType>& other) const;

    BlockVector operator-(const BlockVector<InVectorType>& other) const;

    InVectorType& block(const unsigned int& iblock);

    InVectorType block(const unsigned int& iblock) const;

    BlockVector getSubvector(const unsigned int& ibegin, const unsigned int& iend) const;

    double norm2() const;

    void resize(const unsigned int& nRows);

    void hardCopy(const BlockVector<InVectorType>& other);

    void softCopy(const BlockVector<InVectorType>& other);

    inline unsigned int nRows() const {return M_nRows;}

    void zero();

    InVectorType collapse() const;

protected:

    Grid            M_vectorGrid;
    unsigned int    M_nRows;
};

}

#include "BlockVector_imp.hpp"

#endif // BLOCKVECTOR_HPP
