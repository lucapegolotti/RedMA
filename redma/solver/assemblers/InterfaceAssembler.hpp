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

#ifndef INTERFACEASSEMBLER_HPP
#define INTERFACEASSEMBLER_HPP

#include <redma/RedMA.hpp>

#include <redma/solver/assemblers/aAssembler.hpp>

namespace RedMA
{

template <class InVectorType, class InMatrixType>
class Interface
{
    typedef aAssembler<InVectorType COMMA InMatrixType>         AssemblerType;
public:
    Interface(SHP(AssemblerType) assemblerFather, const unsigned int& indexFather,
              SHP(AssemblerType) assemblerChild, const unsigned int& indexChild);

    SHP(AssemblerType)      M_assemblerFather;
    SHP(AssemblerType)      M_assemblerChild;
    unsigned int            M_indexFather;
    unsigned int            M_indexChild;
};

template <class InVectorType, class InMatrixType>
class InterfaceAssembler
{
public:


protected:

};

}

#include "InterfaceAssembler_imp.hpp"

#endif // INTERFACEASSEMBLER_HPP
