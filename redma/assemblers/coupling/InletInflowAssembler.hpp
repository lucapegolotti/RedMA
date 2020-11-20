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

class InletInflowAssembler  : public InterfaceAssembler
{
    typedef aAssembler         AssemblerType;

public:
    InletInflowAssembler(const DataContainer& data,
                         const Interface& interface);

    virtual void addContributionRhs(const double& time,
                                    shp<BlockVector> rhs,
                                    shp<BlockVector> sol,
                                    const unsigned int& nPrimalBlocks) override;

    virtual void addContributionJacobianRhs(const double& time,
                                            shp<BlockMatrix> jac,
                                            shp<BlockVector> sol,
                                            const unsigned int& nPrimalBlocks) override;

};

}

#endif // INTERFACEASSEMBLER_HPP
