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

#ifndef BCMANAGER_HPP
#define BCMANAGER_HPP

#include <redma/RedMA.hpp>
#include <redma/geometry/TreeStructure.hpp>
#include <redma/solver/array/BlockVector.hpp>
#include <redma/solver/array/BlockMatrix.hpp>
#include <redma/solver/array/VectorEp.hpp>
#include <redma/solver/array/MatrixEp.hpp>

#include <lifev/core/filter/GetPot.hpp>
#include <lifev/core/fem/BCHandler.hpp>

namespace RedMA
{

class BCManager
{
public:
    BCManager(const GetPot& datafile, SHP(TreeNode) treeNode);

    void applyDirichletBCs(const double& time, BlockVector<VectorEp>& input,
                           SHP(FESPACE) fespace, const unsigned int& index) const;

    void apply0DirichletBCs(BlockVector<VectorEp>& input,
                            SHP(FESPACE) fespace,
                            const unsigned int& index) const;

    void apply0DirichletMatrix(BlockMatrix<MatrixEp>& input,
                               SHP(FESPACE) fespace,
                               const unsigned int& index,
                               const double& diagCoefficient) const;

    void setInflow(std::function<double(double)> inflow);

private:
    static double poiseulle(const double& t, const double& x, const double& y,
                            const double& z, const unsigned int& i,
                            const GeometricFace& face,
                            const std::function<double(double)> inflow);

    static double fZero(const double& t, const double& x, const double& y,
                        const double& z, const unsigned int& i);

    SHP(LifeV::BCHandler) createBCHandler0Dirichlet() const;

    SHP(TreeNode)                 M_treeNode;
    GetPot                        M_datafile;
    std::function<double(double)> M_inflow;

    const unsigned int            inletFlag = 1;
    const unsigned int            wallFlag = 10;
    const unsigned int            inletRing = 30;
    const unsigned int            outletRing = 31;
};

}

#endif // BCMANAGER_HPP
