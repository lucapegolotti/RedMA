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
    NavierStokesAssemblerFE(const DataContainer& data, SHP(TreeNode) treeNode,
                            std::string stabilizationName = "");

    void setup() override;

    void addConvectiveMatrixRightHandSide(SHP(aVector) sol,
                                          SHP(aMatrix) mat) override;

    void addConvectiveTermJacobianRightHandSide(SHP(aVector) sol,
                                                SHP(aVector) lifting,
                                                SHP(aMatrix) mat) override;

    SHP(aMatrix) getMass(const double& time,
                         const SHP(aVector)& sol) override;

    SHP(aMatrix) getMassJacobian(const double& time,
                                 const SHP(aVector)& sol) override;

    SHP(aVector) getRightHandSide(const double& time,
                                  const SHP(aVector)& sol) override;

    SHP(aMatrix) getJacobianRightHandSide(const double& time,
                                          const SHP(aVector)& sol) override;


protected:
    SHP(LifeV::Exporter<MESH>)                        M_exporter;
    SHP(VECTOREPETRA)                                 M_velocityExporter;
    SHP(VECTOREPETRA)                                 M_WSSExporter;
    SHP(VECTOREPETRA)                                 M_pressureExporter;
    std::string                                       M_name;
    SHP(BlockVector)                                  M_extrapolatedSolution;
    SHP(NavierStokesStabilization)                    M_stabilization;
    std::string                                       M_stabilizationName;
};

}

#endif // NAVIERSTOKESASSEMBLERFE_HPP
