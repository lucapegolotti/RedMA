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

#ifndef aMATRIX_HPP
#define aMATRIX_HPP

#include <redma/RedMA.hpp>
#include <redma/array/aDataWrapper.hpp>
#include <redma/array/aVector.hpp>
#include <redma/utils/Exception.hpp>
#include <redma/array/TypesUtils.hpp>

#include <memory>

namespace RedMA
{

class aMatrix : public aDataWrapper
{
public:

    aMatrix() : M_nRows(0), M_nCols(0) {}

    virtual ~aMatrix() {};

    virtual void add(shp<aMatrix> other) = 0;

    virtual void multiplyByScalar(const double& coeff) = 0;

    virtual shp<aVector> multiplyByVector(shp<aVector> vector) = 0;

    virtual shp<aMatrix> multiplyByMatrix(shp<aMatrix> other) = 0;

    virtual shp<aMatrix> transpose() const = 0;

    virtual void dump(std::string filename) const = 0;

    inline unsigned int nRows() const {return M_nRows;}

    inline unsigned int nCols() const {return M_nCols;}

    virtual shp<aMatrix> block(const unsigned int& row, const unsigned int& col) const
    {
        throw new Exception("block(row,col) function not overloaded!");
    }

protected:

    unsigned int        M_nRows;
    unsigned int        M_nCols;
};

}

#endif // aMATRIX_HPP
