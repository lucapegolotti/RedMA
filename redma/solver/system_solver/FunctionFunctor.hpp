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

#ifndef FUNCTIONFUNCTOR_HPP
#define FUNCTIONFUNCTOR_HPP

#include <redma/RedMA.hpp>
#include <redma/array/BlockVector.hpp>

#include <functional>

namespace RedMA
{

template <class InputType, class OutputType>
class FunctionFunctor
{
public:
    FunctionFunctor(std::function<OutputType(InputType)> fct) : M_function(fct) {}

    inline OutputType operator()(const InputType& input) const
    {
        return M_function(input);
    }

private:
    std::function<OutputType(InputType)> M_function;
};

}

#endif // FUNCTIONFUNCTOR_HPP
