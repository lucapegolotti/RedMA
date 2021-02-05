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
#include <redma/array/BlockVector.hpp>
#include <redma/array/BlockMatrix.hpp>
#include <redma/array/DistributedVector.hpp>
#include <redma/array/SparseMatrix.hpp>
#include <redma/problem/DataContainer.hpp>
#include <redma/solver/time_marching_algorithms/TimeMarchingAlgorithmFactory.hpp>
#include <redma/boundary_conditions/WindkesselModel.hpp>

#include <lifev/core/fem/BCHandler.hpp>

namespace RedMA
{

class BCManager
{
    typedef aTimeMarchingAlgorithm  TimeMarchingAlgorithm;
public:
    BCManager(const DataContainer& datafile, shp<TreeNode> treeNode);

    void applyDirichletBCs(const double& time, BlockVector& input,
                           shp<FESPACE> fespace, const unsigned int& index,
                           const bool& ringOnly = false) const;

    void apply0DirichletBCs(BlockVector& input,
                            shp<FESPACE> fespace,
                            const unsigned int& index,
                            const bool& ringOnly = false) const;

    void apply0DirichletMatrix(BlockMatrix& input,
                               shp<FESPACE> fespace,
                               const unsigned int& index,
                               const double& diagCoefficient,
                               const bool& ringOnly = false) const;

    void applyInflowDirichletBCs(shp<LifeV::BCHandler> bcs, const bool& zeroFlag = false) const;

    void applyInflowNeumannBCs(shp<LifeV::BCHandler> bcs, const bool& zeroFlag = false) const;

    double getOutflowNeumannBC(const double& time, const double& flag, const double& rate);

    // actually derivative wrt to flowrate
    double getOutflowNeumannJacobian(const double& time, const double& flag, const double& rate);

    void postProcess();

    shp<VECTOREPETRA> computeBoundaryIndicator(shp<FESPACE> fespace);

    inline bool useStrongDirichlet() const {return M_strongDirichlet;}

private:
    static double poiseuilleInflow(const double& t, const double& x, const double& y,
                                  const double& z, const unsigned int& i,
                                  const GeometricFace& face,
                                  const std::function<double(double)> inflow,
                                  const double& coefficient);

    static double neumannInflow(const double& t, const double& x, const double& y,
                                const double& z, const unsigned int& i,
                                std::function<double(double)> inflowLaw);

    void checkInflowLaw();

    static double fZero(const double& t, const double& x, const double& y,
                        const double& z, const unsigned int& i);

    static double fOne(const double& t, const double& x, const double& y,
                       const double& z, const unsigned int& i);

    static double constantFunction(const double& t, const double& x, const double& y,
                                   const double& z, const unsigned int& i,
                                   const double& K);

    shp<LifeV::BCHandler> createBCHandler0Dirichlet(const bool& ringOnly = false) const;

    void addInletBC(shp<LifeV::BCHandler> bcs,
                    const bool& ringOnly = false,
                    const bool& zeroFlag = false) const;

    void parseOutflowNeumannData();

    shp<TreeNode>                                    M_treeNode;
    DataContainer                                    M_data;

    std::string                                      M_inletBCType;
    std::function<double(double)>                    M_inflow;
    bool                                             M_strongDirichlet;
    double                                           M_coefficientInflow;

    unsigned int                                     M_inletFlag;
    unsigned int                                     M_wallFlag;
    unsigned int                                     M_inletRingFlag;
    unsigned int                                     M_outletRingFlag;

    // key is the outlet index (more than one for bifurcations)
    std::map<unsigned int, shp<WindkesselModel>>     M_models;
};

}

#endif // BCMANAGER_HPP
