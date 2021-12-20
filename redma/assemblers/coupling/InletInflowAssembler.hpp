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

#ifndef INLETINFLOWASSEMBLER_HPP
#define INLETINFLOWASSEMBLER_HPP

#include <redma/assemblers/coupling/InterfaceAssembler.hpp>

namespace RedMA
{

/// \brief Class for the assembly of the matrix handling the weak dirichlet bcs.
class InletInflowAssembler  : public InterfaceAssembler
{
    typedef aAssembler         AssemblerType;

public:
    /*! \brief Constructor taking a DataContainer and an Interface as arguments.
     *
     * \param data The DataContainer of the problem.
     * \param interface The interface.
     * \param addNoSlipBC True in no-slip BCs at the vessel wall are desired (default)
     */
    InletInflowAssembler(const DataContainer& data,
                         const Interface& interface,
                         const bool& addNoSlipBC = true);

    /*! \brief Add coupling contribution to a right-hand side.
    *
    * \param time Current time (needed in derived classes).
    * \param rhs Shared pointer to the right-hand side.
    * \param sol Shared pointer to the solution.
    * \param nPrimalBlocks number of primal blocks in the problem.
    */
    virtual void addContributionRhs(const double& time,
                                    shp<BlockVector> rhs,
                                    shp<BlockVector> sol,
                                    const unsigned int& nPrimalBlocks) override;

    /*! \brief Add coupling contribution to the jacobian right-hand side.
    *
    * \param time Current time (needed in derived classes).
    * \param jac Shared pointer to the jacobian of the right-hand side.
    * \param sol Shared pointer to the solution.
    * \param nPrimalBlocks Number of primal blocks in the problem.
    */
    virtual void addContributionJacobianRhs(const double& time,
                                            shp<BlockMatrix> jac,
                                            shp<BlockVector> sol,
                                            const unsigned int& nPrimalBlocks) override;

};

}

#endif // INLETINFLOWASSEMBLER_HPP
