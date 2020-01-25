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
#include <redma/geometry/TreeStructure.hpp>

namespace RedMA
{

template <class InVectorType, class InMatrixType>
class BlockAssembler : public aAssembler<InVectorType, InMatrixType>
{
public:
    BlockAssembler(const GetPot& datafile, SHP(TreeStructure) tree);

    virtual void exportSolution(const double& t);

    virtual void postProcess();

    virtual BlockMatrix<InMatrixType> getMass(const double& time,
                                      const BlockVector<InVectorType>& sol);

    virtual BlockVector<InVectorType> getRightHandSide(const double& time,
                                      const BlockVector<InVectorType>& sol);

    virtual BlockMatrix<InMatrixType> getJacobianRightHandSide(const double& time,
                                      const BlockVector<InVectorType>& sol);

protected:
    GetPot              M_datafile;
    SHP(TreeStructure)  M_tree;
};

}

#include "BlockAssembler_imp.hpp"

#endif // BLOCKASSEMBLER_HPP
