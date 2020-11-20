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

#ifndef NAVIERSTOKESMODEL_HPP
#define NAVIERSTOKESMODEL_HPP

#include <redma/assemblers/models/StokesModel.hpp>
// #include <redma/solver/finite_element/SUPGStabilization.hpp>
// #include <redma/solver/assemblers/VMSStabilization.hpp>
// #include <redma/solver/assemblers/HFStabilization.hpp>

namespace RedMA
{

class NavierStokesModel
{
public:
    NavierStokesModel(const DataContainer& data, shp<TreeNode> treeNode);

    // virtual void setup() override;

    // virtual BlockMatrix getMass(const double& time,
    //                             const BlockVector& sol) override;
    //
    // virtual BlockMatrixgetMassJacobian(const double& time,
    //                                    const BlockVector& sol) override;
    //
    // virtual BlockVector getRightHandSide(const double& time,
    //                                      const BlockVector& sol) override;
    //
    // virtual BlockMatrix getJacobianRightHandSide(const double& time,
    //                                              const BlockVector& sol) override;
    //
    // virtual BlockVector getNonLinearTerm() override {return M_nonLinearTerm;}
    //
    // virtual void RBsetup() override;

    virtual void addConvectiveMatrixRightHandSide(shp<aVector> sol,
                                                  shp<aMatrix> mat) = 0;

    virtual void addConvectiveTermJacobianRightHandSide(shp<aVector> sol,
                                                        shp<aVector> lifting,
                                                        shp<aMatrix> mat) = 0;

protected:

    // shp<BlockVector)                            M_nonLinearTerm;
    //
    // // this is relative to the rb part
    // std::vector<std::vector<BlockVector>>       M_nonLinearTermsDecomposition;
};

}

#endif // NAVIERSTOKESMODEL_HPP
