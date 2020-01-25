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

#ifndef aASSEMBLER_HPP
#define aASSEMBLER_HPP

#include <redma/RedMA.hpp>
#include <redma/solver/array/BlockMatrix.hpp>
#include <redma/geometry/TreeStructure.hpp>

#include <lifev/core/filter/GetPot.hpp>

namespace RedMA
{

template <class InVectorType, class InMatrixType>
class aAssembler
{
public:
    aAssembler(const GetPot& datafile);

    aAssembler(const GetPot& datafile, SHP(TreeNode) node);

    virtual void setup() = 0;

    virtual void exportSolution(const double& t) = 0;

    virtual void postProcess() = 0;

    virtual BlockMatrix<InMatrixType> getMass(const double& time,
                                      const BlockVector<InVectorType>& sol) = 0;

    virtual BlockVector<InVectorType> getRightHandSide(const double& time,
                                      const BlockVector<InVectorType>& sol) = 0;

    virtual BlockMatrix<InMatrixType> getJacobianRightHandSide(const double& time,
                                      const BlockVector<InVectorType>& sol) = 0;

protected:
    GetPot                M_datafile;
    SHP(TreeNode)         M_treeNode;
};

}

#include "aAssembler_imp.hpp"

#endif // aASSEMBLER_HPP
