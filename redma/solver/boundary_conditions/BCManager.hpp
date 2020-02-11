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
#include <redma/solver/problem/DataContainer.hpp>
#include <redma/solver/time_marching_algorithms/TimeMarchingAlgorithmFactory.hpp>
#include <redma/solver/boundary_conditions/WindkesselModel.hpp>

#include <lifev/core/fem/BCHandler.hpp>

namespace RedMA
{

class BCManager
{
    typedef aTimeMarchingAlgorithm<Double COMMA Double>  TimeMarchingAlgorithm;
public:
    BCManager(const DataContainer& datafile, SHP(TreeNode) treeNode);

    void applyDirichletBCs(const double& time, BlockVector<VectorEp>& input,
                           SHP(FESPACE) fespace, const unsigned int& index) const;

    void apply0DirichletBCs(BlockVector<VectorEp>& input,
                            SHP(FESPACE) fespace,
                            const unsigned int& index) const;

    void apply0DirichletMatrix(BlockMatrix<MatrixEp>& input,
                               SHP(FESPACE) fespace,
                               const unsigned int& index,
                               const double& diagCoefficient) const;

    double getNeumannBc(const double& time, const double& flag, const double& rate);

    // actually derivative wrt to flowrate
    double getNeumannJacobian(const double& time, const double& flag, const double& rate);

    void postProcess();

    inline bool useStrongDirichlet() const {return M_strongDirichlet;}

private:
    static double poiseuille(const double& t, const double& x, const double& y,
                             const double& z, const unsigned int& i,
                             const GeometricFace& face,
                             const std::function<double(double)> inflow,
                             const double& coefficient);

    static double fZero(const double& t, const double& x, const double& y,
                        const double& z, const unsigned int& i);

    static double constantFunction(const double& t, const double& x, const double& y,
                                   const double& z, const unsigned int& i,
                                   const double& K);

    static double fZero2(double t);

    SHP(LifeV::BCHandler) createBCHandler0Dirichlet() const;

    void addInletBC(SHP(LifeV::BCHandler) bcs,
                    std::function<double(double)> law) const;

    void parseNeumannData();

    SHP(TreeNode)                                    M_treeNode;
    DataContainer                                    M_data;
    std::function<double(double)>                    M_inflow;
    std::function<double(double)>                    M_inflowDt;
    bool                                             M_strongDirichlet;

    const unsigned int                               inletFlag = 1;
    const unsigned int                               wallFlag = 10;
    const unsigned int                               inletRing = 30;
    const unsigned int                               outletRing = 31;

    // key is the outlet index (more than one for bifurcations)
    std::map<unsigned int,SHP(WindkesselModel)>      M_models;

    double                                           M_coefficientInflow;
};

}

#endif // BCMANAGER_HPP
