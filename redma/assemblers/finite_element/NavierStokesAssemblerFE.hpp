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

#ifndef NAVIERSTOKESASSEMBLERFE_HPP
#define NAVIERSTOKESASSEMBLERFE_HPP

#include <redma/assemblers/finite_element/StokesAssemblerFE.hpp>
#include <redma/assemblers/models/NavierStokesModel.hpp>
#include <redma/assemblers/finite_element/SUPGStabilization.hpp>

namespace RedMA
{

class NavierStokesAssemblerFE : public StokesAssemblerFE, public NavierStokesModel
{
public:
    NavierStokesAssemblerFE(const DataContainer& data, shp<TreeNode> treeNode,
                            std::string stabilizationName = "");

    void setup() override;

    void addConvectiveMatrixRightHandSide(shp<aVector> sol,
                                          shp<aMatrix> mat) override;

    void addConvectiveTermJacobianRightHandSide(shp<aVector> sol,
                                                shp<aVector> lifting,
                                                shp<aMatrix> mat) override;

    shp<aMatrix> getMass(const double& time,
                         const shp<aVector>& sol) override;

    shp<aMatrix> getMassJacobian(const double& time,
                                 const shp<aVector>& sol) override;

    shp<aVector> getRightHandSide(const double& time,
                                  const shp<aVector>& sol) override;

    shp<aMatrix> getJacobianRightHandSide(const double& time,
                                          const shp<aVector>& sol) override;


protected:
    shp<NavierStokesStabilization>                    M_stabilization;
    std::string                                       M_stabilizationName;
};

}

#endif // NAVIERSTOKESASSEMBLERFE_HPP
