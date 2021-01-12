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

#ifndef aASSEMBLERRB_HPP
#define aASSEMBLERRB_HPP

#include <redma/assemblers/abstract/aAssembler.hpp>

namespace RedMA
{

/// Abstract assembler class for a reduced basis (RB) problem.
class aAssemblerRB : public aAssembler
{
public:
    /*! \brief Constructor taking a datafile as argument.
     *
     * \param datafile The datafile.
     */
    aAssemblerRB(const DataContainer& datafile) :
      aAssembler(datafile)
    {}

    /*! \brief Constructor taking a datafile and a TreeNode as argument.
     *
     * \param datafile The datafile.
     * \param datafile The TreeNode encoding the physical domain.
     */
    aAssemblerRB(const DataContainer& datafile, shp<TreeNode> node) :
      aAssembler(datafile, node)
    {}

    /*! \brief Project the RB function onto the finite element space.
     *
     * \param rbSolution RB coefficients of the solutions
     * \return Shared pointer to aVector of the projected function.
     */
    virtual shp<aVector> convertFunctionRBtoFEM(shp<aVector> rbSolution) const = 0;
};

}

#endif // aASSEMBLERRB_HPP
