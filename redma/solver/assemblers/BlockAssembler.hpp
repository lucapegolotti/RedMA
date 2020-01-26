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

#ifndef BLOCKASSEMBLER_HPP
#define BLOCKASSEMBLER_HPP

#include <redma/RedMA.hpp>
#include <redma/solver/assemblers/aAssembler.hpp>
#include <redma/solver/assemblers/AssemblerFactory.hpp>
#include <redma/solver/assemblers/InterfaceAssembler.hpp>

#include <redma/geometry/TreeStructure.hpp>

namespace RedMA
{

template <class InVectorType, class InMatrixType>
class BlockAssembler : public aAssembler<InVectorType, InMatrixType>
{
    typedef typename InVectorType::InnerType              VInner;
    typedef typename InMatrixType::InnerType              MInner;

public:
    BlockAssembler(const GetPot& datafile, const TreeStructure& tree);

    virtual void setup() override;

    virtual void exportSolution(const double& t,
                                const BlockVector<InVectorType>& sol) override;

    virtual void postProcess() override;

    virtual BlockMatrix<InMatrixType> getMass(const double& time,
                                      const BlockVector<InVectorType>& sol) override;

    virtual BlockVector<InVectorType> getRightHandSide(const double& time,
                                      const BlockVector<InVectorType>& sol) override;

    virtual BlockMatrix<InMatrixType> getJacobianRightHandSide(const double& time,
                                      const BlockVector<InVectorType>& sol) override;

protected:
    GetPot                                                        M_datafile;
    TreeStructure                                                 M_tree;
    std::map<unsigned int, SHP(aAssembler<VInner COMMA MInner>)>  M_primalAssemblers;
    std::vector<SHP(InterfaceAssembler<VInner COMMA MInner>)>     M_dualAssemblers;
    unsigned int                                                  M_numberBlocks;
};

}

#include "BlockAssembler_imp.hpp"

#endif // BLOCKASSEMBLER_HPP
