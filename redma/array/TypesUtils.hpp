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

#ifndef TYPESUTILS_HPP
#define TYPESUTILS_HPP

#include <redma/utils/Exception.hpp>
#include <memory>

#define ZEROTHRESHOLD       1e-100

namespace RedMA
{

/// Convert a class to another type, if possible.
template <class DType, class ContainerType>
shp<DType> convert(shp<ContainerType> container)
{
    if (container->type() != DType().type())
    {
        std::string msg = "Error in convert: converting ";
        msg += std::to_string(container->type());
        msg += " to ";
        msg += std::to_string(DType().type());
        msg += "\n";
        throw new Exception(msg);
    }
    return spcast<DType>(container);
}

}

#endif // TYPESUTILS_HPP
