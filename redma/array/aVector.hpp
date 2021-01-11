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

#ifndef aVECTOR_HPP
#define aVECTOR_HPP

#include <redma/RedMA.hpp>
#include <redma/array/aDataWrapper.hpp>
#include <redma/array/Datatypes.hpp>
#include <redma/utils/Exception.hpp>
#include <redma/array/TypesUtils.hpp>

#include <memory>

namespace RedMA
{

class aVector : public aDataWrapper
{
public:
    aVector() : M_nRows(0) {}

    virtual ~aVector() {};

    virtual void add(shp<aVector> other) = 0;

    virtual void multiplyByScalar(const double& coeff) = 0;

    virtual std::string getString(const char& delimiter) const = 0;

    virtual double norm2() const = 0;

    virtual double operator()(unsigned int index) = 0;

    inline unsigned int nRows() const {return M_nRows;}

protected:
    unsigned int        M_nRows;
};

}

#endif // aVECTOR_HPP
