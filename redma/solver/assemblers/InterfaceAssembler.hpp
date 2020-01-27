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
#include <redma/solver/array/BlockMatrix.hpp>
#include <redma/solver/array/BlockVector.hpp>
#include <redma/solver/array/MatrixEp.hpp>
#include <redma/solver/array/VectorEp.hpp>

namespace RedMA
{

template <class InVectorType, class InMatrixType>
class Interface
{
    typedef aAssembler<InVectorType COMMA InMatrixType>         AssemblerType;
public:
    Interface(SHP(AssemblerType) assemblerFather, const unsigned int& indexFather,
              SHP(AssemblerType) assemblerChild, const unsigned int& indexChild,
              const unsigned int& interfaceID);

    SHP(AssemblerType)      M_assemblerFather;
    SHP(AssemblerType)      M_assemblerChild;
    unsigned int            M_indexFather;
    unsigned int            M_indexChild;
    unsigned int            M_ID;
};

template <class InVectorType, class InMatrixType>
class InterfaceAssembler
{
public:
    InterfaceAssembler(const Interface<InVectorType, InMatrixType>& interface);

    void setup();

    void buildCouplingMatrices();

    void addContributionRhs(BlockVector<BlockVector<InVectorType>>& rhs,
                            const BlockVector<BlockVector<InVectorType>>& sol,
                            const unsigned int& nPrimalBlocks);

    void addContributionJacobianRhs(BlockMatrix<BlockMatrix<InMatrixType>>& jac,
                                    const BlockVector<BlockVector<InVectorType>>& sol,
                                    const unsigned int& nPrimalBlocks);

    inline Interface<InVectorType, InMatrixType> getInterface() const {return M_interface;};

private:
    Interface<InVectorType, InMatrixType>           M_interface;

    BlockMatrix<InMatrixType>                       M_fatherBT;
    BlockMatrix<InMatrixType>                       M_fatherB;
    BlockMatrix<InMatrixType>                       M_childBT;
    BlockMatrix<InMatrixType>                       M_childB;
};

}

#include "InterfaceAssembler_imp.hpp"

#endif // INTERFACEASSEMBLER_HPP