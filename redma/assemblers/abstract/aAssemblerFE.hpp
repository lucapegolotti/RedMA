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

#ifndef aASSEMBLERFE_HPP
#define aASSEMBLERFE_HPP

#include <redma/assemblers/abstract/aAssembler.hpp>

namespace RedMA
{

class aAssemblerFE : public aAssembler
{
public:
    virtual SparseMatrix getNorm(const unsigned int& fieldIndex, bool bcs = true) {return SparseMatrix();}

    virtual SparseMatrix getConstraintMatrix() {return SparseMatrix();}

};

}

#endif // aASSEMBLERFE_HPP