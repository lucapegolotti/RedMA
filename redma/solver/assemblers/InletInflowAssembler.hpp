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

#include <redma/solver/assemblers/InterfaceAssembler.hpp>

namespace RedMA
{

template <class InVectorType, class InMatrixType>
class InletInflowAssembler : public InterfaceAssembler<InVectorType, InMatrixType>
{
    typedef aAssembler<InVectorType COMMA InMatrixType>         AssemblerType;

public:
    InletInflowAssembler(const DataContainer& data,
                         const Interface<InVectorType, InMatrixType>& interface);

    virtual void addContributionRhs(const double& time,
                                    BlockVector<BlockVector<InVectorType>>& rhs,
                                    const BlockVector<BlockVector<InVectorType>>& sol,
                                    const unsigned int& nPrimalBlocks);

    virtual void addContributionJacobianRhs(const double& time,
                                            BlockMatrix<BlockMatrix<InMatrixType>>& jac,
                                            const BlockVector<BlockVector<InVectorType>>& sol,
                                            const unsigned int& nPrimalBlocks);

};

}

#include "InletInflowAssembler_imp.hpp"

#endif // INTERFACEASSEMBLER_HPP
