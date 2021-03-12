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

#ifndef REDMA_MEMBRANEASSEMBLERFE_HPP
#define REDMA_MEMBRANEASSEMBLERFE_HPP

#include <redma/RedMA.hpp>
#include <redma/assemblers/finite_element/NavierStokesAssemblerFE.hpp>
#include <redma/solver/time_marching_algorithms/TimeMarchingAlgorithmFactory.hpp>
#include <redma/geometry/MembraneThicknessComputer.hpp>

namespace RedMA
{

class MembraneAssemblerFE : public NavierStokesAssemblerFE {

public:

    MembraneAssemblerFE(const DataContainer& datafile,
                        shp<TreeNode> treeNode,
                        std::string stabilizationName = "");

    void setup() override;

    shp<aMatrix> assembleMass(shp<BCManager> bcManager) override;

    shp<aMatrix> assembleBoundaryMass(shp<BCManager> bcManager, bool verbose = false);

    shp<aMatrix> assembleWallBoundaryMass(shp<BCManager> bcManager, bool verbose = false);

    shp<aMatrix> assembleStiffness(shp<BCManager> bcManager) override;

    shp<aMatrix> assembleBoundaryStiffness(shp<BCManager> bcManager, bool verbose = false);

    shp<aVector> getRightHandSide(const double& time,
                                  const shp<aVector>& sol) override;

    shp<aMatrix> getJacobianRightHandSide(const double& time,
                                          const shp<aVector>& sol) override;

    std::vector<shp<aMatrix>> getBoundaryMatrices() const;

    void postProcess(const double& t, const shp<aVector>& sol) override;

    void setExporter() override;

    void setDisplacementExporter(shp<VECTOREPETRA> displacement);

    inline double getWallElasticity() const {return M_wall_elasticity;}

    inline double getWallViscoelasticity() const {return M_wall_viscoelasticity;}

    inline shp<VECTOREPETRA> getBoundaryIndicator() const {return M_boundaryIndicator;}

protected:

    void computeLameConstants();

    void computeThickness();

    void exportThickness();

    shp<aTimeMarchingAlgorithm>                     M_TMA_Displacements;

    double                                          M_lameI;
    double                                          M_lameII;
    double                                          M_membrane_density;
    double                                          M_transverse_shear_coeff;
    double                                          M_wall_elasticity;
    double                                          M_wall_viscoelasticity;
    unsigned int                                    M_wallFlag;

    shp<BlockMatrix>                                M_boundaryStiffness;
    shp<BlockMatrix>                                M_boundaryMass;
    shp<BlockMatrix>                                M_wallBoundaryMass;

    shp<VECTOREPETRA>                               M_displacementExporter;
    shp<VECTOREPETRA>                               M_boundaryIndicator;

};

};


#endif //REDMA_MEMBRANEASSEMBLERFE_HPP
