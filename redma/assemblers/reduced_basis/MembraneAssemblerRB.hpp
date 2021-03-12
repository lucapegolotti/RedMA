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

#ifndef REDMA_MEMBRANEASSEMBLERRB_HPP
#define REDMA_MEMBRANEASSEMBLERRB_HPP

#include <redma/assemblers/reduced_basis/NavierStokesAssemblerRB.hpp>
#include <redma/assemblers/finite_element/MembraneAssemblerFE.hpp>
#include <redma/solver/time_marching_algorithms/TimeMarchingAlgorithmFactory.hpp>

namespace RedMA {

    class MembraneAssemblerRB : public NavierStokesAssemblerRB {

    public:
        MembraneAssemblerRB(const DataContainer &data, shp<TreeNode> treeNode);

        void RBsetup() override;

        shp<aVector> getRightHandSide(const double &time,
                                      const shp<aVector> &sol) override;

        shp<aMatrix> getJacobianRightHandSide(const double& time,
                                              const shp<aVector>& sol) override;

        void postProcess(const double& t, const shp<aVector>& sol) override;

    protected:

        inline shp<MembraneAssemblerFE> getFEAssembler()
        {
            return getFEAssemblerAs<MembraneAssemblerFE>();
        }

        shp<BlockMatrix>                                M_reducedBoundaryMass;
        shp<BlockMatrix>                                M_reducedBoundaryStiffness;
        shp<BlockMatrix>                                M_reducedWallBoundaryMass;

        shp<aTimeMarchingAlgorithm>                     M_TMA_Displacements;
    };

} // Namespace RedMA

#endif //REDMA_MEMBRANEASSEMBLERRB_HPP
