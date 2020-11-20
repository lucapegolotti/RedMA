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

#ifndef NAVIERSTOKESASSEMBLERRB_HPP
#define NAVIERSTOKESASSEMBLERRB_HPP

#include <redma/assemblers/abstract/aAssemblerFE.hpp>
#include <redma/assemblers/reduced_basis/StokesAssemblerRB.hpp>
#include <redma/assemblers/models/NavierStokesModel.hpp>
#include <redma/assemblers/abstract/aAssemblerRB.hpp>
#include <redma/reduced_basis/RBBases.hpp>

namespace RedMA
{

class NavierStokesAssemblerRB : public StokesAssemblerRB, public NavierStokesModel
{
public:
    NavierStokesAssemblerRB(const DataContainer& data, shp<TreeNode> treeNode);

    void addConvectiveMatrixRightHandSide(shp<aVector> sol,
                                          shp<aMatrix> mat) override;

    void addConvectiveTermJacobianRightHandSide(shp<aVector> sol,
                                                shp<aVector> lifting,
                                                shp<aMatrix> mat) override;

    shp<aVector> getRightHandSide(const double& time,
                                  const shp<aVector>& sol) override;

    shp<aMatrix> getJacobianRightHandSide(const double& time,
                                          const shp<aVector>& sol) override;

    virtual void RBsetup() override;

protected:
    std::vector<std::vector<shp<BlockVector>>>       M_nonLinearTermsDecomposition;
    shp<BlockVector>                                 M_nonLinearTerm;
    bool                                             M_exactJacobian;
};

}

#endif // NAVIERSTOKESASSEMBLERRB_HPP
