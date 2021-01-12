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

/// Abstract assembler class for a finite element (FE) problem.
class aAssemblerFE : public aAssembler
{
public:
    /*! \brief Constructor taking a datafile as argument.
     *
     * \param datafile The datafile.
     */
    aAssemblerFE(const DataContainer& datafile) :
      aAssembler(datafile)
    {}

    /*! \brief Constructor taking a datafile and a TreeNode as argument.
     *
     * \param datafile The datafile.
     * \param datafile The TreeNode encoding the physical domain.
     */
    aAssemblerFE(const DataContainer& datafile, shp<TreeNode> node) :
      aAssembler(datafile, node)
    {}

    /*! \brief Getter for the norm corresponding to a specific field.
     *
     * \param fieldIndex Index of the desired field.
     * \param bcs If true, applies bcs to the norm matrix.
     * \return Shared pointer to aMatrix of the desired norm.
     */
    virtual shp<aMatrix> getNorm(const unsigned int& fieldIndex, bool bcs = true) {return shp<SparseMatrix>();}

    /*! \brief Getter for the constraint matrix (e.g., divergence matrix in Stokes).
     *
     * If not overloaded, returns an empty SparseMatrix.
     *
     * \return Shared pointer to aMatrix of the constraint matrix.
     */
    virtual shp<aMatrix> getConstraintMatrix() {return shp<SparseMatrix>();}

    /*! \brief Getter for the FE lifting.
     *
     * \param time Current time.
     * \return Shared pointer to aVector of the FE lifting.
     */
    virtual shp<aVector> getFELifting(const double& time) const = 0;
};

}

#endif // aASSEMBLERFE_HPP
