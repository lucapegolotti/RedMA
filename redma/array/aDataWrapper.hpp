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

#ifndef aDATAWRAPPER_HPP
#define aDATAWRAPPER_HPP

#include <redma/RedMA.hpp>
#include <redma/array/Datatypes.hpp>
#include <redma/utils/Exception.hpp>
#include <redma/array/TypesUtils.hpp>

#include <memory>

namespace RedMA
{

class aDataWrapper
{
public:

    virtual ~aDataWrapper() {}

    virtual Datatype type() const = 0;

    virtual shp<void> data() const = 0;

    virtual void setData(shp<void> data) = 0;

    virtual bool isZero() const = 0;

    virtual aDataWrapper* clone() const = 0;

    virtual void shallowCopy(shp<aDataWrapper> other) = 0;

    virtual void deepCopy(shp<aDataWrapper> other) = 0;

    virtual void dump(std::string namefile) const = 0;
};

}

#endif // aDATAWRAPPER_HPP
