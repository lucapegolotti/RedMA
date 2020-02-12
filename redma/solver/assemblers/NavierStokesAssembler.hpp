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

#ifndef NAVIERSTOKESASSEMBLER_HPP
#define NAVIERSTOKESASSEMBLER_HPP

#include <redma/solver/assemblers/StokesAssembler.hpp>
#include <redma/solver/assemblers/SUPGStabilization.hpp>
#include <redma/solver/assemblers/VMSStabilization.hpp>

namespace RedMA
{

template <class InVectorType, class InMatrixType>
class NavierStokesAssembler : public StokesAssembler<InVectorType, InMatrixType>
{
public:
    NavierStokesAssembler(const DataContainer& data, SHP(TreeNode) treeNode);

    virtual void setup() override;

    virtual BlockMatrix<InMatrixType> getMass(const double& time,
                                      const BlockVector<InVectorType>& sol) override;

    virtual BlockMatrix<InMatrixType> getMassJacobian(const double& time,
                                      const BlockVector<InVectorType>& sol) override;

    virtual BlockVector<InVectorType> getRightHandSide(const double& time,
                                      const BlockVector<InVectorType>& sol) override;

    virtual BlockMatrix<InMatrixType> getJacobianRightHandSide(const double& time,
                                      const BlockVector<InVectorType>& sol) override;

protected:
    void addConvectiveMatrixRightHandSide(const BlockVector<InVectorType>& sol,
                                          BlockMatrix<InMatrixType>& mat);

    void addConvectiveTermJacobianRightHandSide(const BlockVector<InVectorType>& sol,
                                                const BlockVector<InVectorType>& lifting,
                                                BlockMatrix<InMatrixType>& rhs);

    bool                                                     M_useStabilization;
    SHP(NavierStokesStabilization)                           M_stabilization;
};

}

#include "NavierStokesAssembler_imp.hpp"

#endif // NAVIERSTOKESASSEMBLER_HPP
