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
#include <lifev/core/array/VectorSmall.hpp>
#include <lifev/core/array/MatrixSmall.hpp>

namespace RedMA
{

class BCManager
{
    typedef aTimeMarchingAlgorithm    TimeMarchingAlgorithm;

    typedef LifeV::VectorSmall<3>     Vector3D;
    typedef LifeV::MatrixSmall<3,3>   Matrix3D;

public:
    BCManager(const DataContainer& datafile, shp<TreeNode> treeNode);

    void applyDirichletBCs(const double& time, BlockVector& input,
                           shp<FESPACE> fespace, const unsigned int& index,
                           const bool& ringOnly = false);

    void apply0DirichletBCs(BlockVector& input,
                            shp<FESPACE> fespace,
                            const unsigned int& index,
                            const bool& ringOnly = false);

    void apply0DirichletMatrix(BlockMatrix& input,
                               shp<FESPACE> fespace,
                               const unsigned int& index,
                               const double& diagCoefficient,
                               const bool& ringOnly = false);

    void applyInflowDirichletBCs(shp<LifeV::BCHandler> bcs, const bool& zeroFlag = false) const;

    void applyInflowNeumannBCs(shp<LifeV::BCHandler> bcs, const bool& zeroFlag = false) const;

    void applyOutflowNeumannBCs(shp<LifeV::BCHandler> bcs, const bool& zeroFlag = true) const;

    double getOutflowNeumannBC(const double& time, const double& flag, const double& rate);

    // actually derivative wrt to flowrate
    double getOutflowNeumannJacobian(const double& time, const double& flag, const double& rate);

    void postProcess();

    inline bool useStrongDirichlet() const {return M_strongDirichlet;}

    inline std::string getInletBCType() const {return M_inletBCType;}

    std::vector<unsigned int> getWallFlags(const bool& withRings = true) const;

    shp<VECTOREPETRA> computeBoundaryIndicator(shp<FESPACE> fespace, const std::vector<unsigned int> flags) const;

    static double fZero(const double& t, const double& x, const double& y,
                        const double& z, const unsigned int& i);

    static double fOne(const double& t, const double& x, const double& y,
                       const double& z, const unsigned int& i);

    static double constantFunction(const double& t, const double& x, const double& y,
                                   const double& z, const unsigned int& i,
                                   const double& K);

private:
    static double poiseuilleInflow(const double& t, const double& x, const double& y,
                                   const double& z, const unsigned int& i,
                                   const GeometricFace& face,
                                   const std::function<double(double)> inflow,
                                   const double& coefficient);

    static double neumannInflow(const double& t, const double& x, const double& y,
                                const double& z, const unsigned int& i,
                                std::function<double(double)> inflowLaw);

    Matrix3D computeRotationMatrix(Vector3D vec, const unsigned int& index = 2) const;

    std::map<unsigned int, Matrix3D> computeRotationMatrices() const;

    void computeGlobalRotationMatrix(shp<FESPACE> fespace);

    void shiftToNormalTangentialCoordSystem(shp<MATRIXEPETRA> mat, shp<VECTOREPETRA> vec,
                                            shp<FESPACE> fespace);

    void shiftToCartesianCoordSystem(shp<MATRIXEPETRA> mat, shp<VECTOREPETRA> vec,
                                     shp<FESPACE> fespace);

    shp<VECTOREPETRA> computeRingsIndicator(shp<FESPACE> fespace) const;

    std::pair<Vector3D, Vector3D> computeTangentVersors(const Vector3D& normal) const;

    void checkInflowLaw();

    shp<LifeV::BCHandler> createBCHandler0Dirichlet(const bool& ringOnly = false) const;

    shp<LifeV::BCHandler> createBCHandler0DirichletRing() const;

    void addInletBC(shp<LifeV::BCHandler> bcs,
                    const bool& ringOnly = false,
                    const bool& zeroFlag = false) const;

    void parseOutflowNeumannData();

    shp<TreeNode>                                    M_treeNode;
    DataContainer                                    M_data;

    std::string                                      M_inletBCType;
    std::string                                      M_ringConstraint;
    std::function<double(double)>                    M_inflow;
    bool                                             M_strongDirichlet;
    double                                           M_coefficientInflow;

    unsigned int                                     M_inletFlag;
    std::vector<unsigned int>                        M_outletFlags;
    std::vector<unsigned int>                        M_trueOutletFlags;
    unsigned int                                     M_wallFlag;
    unsigned int                                     M_inletRingFlag;
    std::vector<unsigned int>                        M_outletRingFlags;
    std::vector<unsigned int>                        M_trueOutletRingFlags;

    shp<MATRIXEPETRA>                                M_globalRotationMatrix;

    // key is the outlet index (more than one for bifurcations)
    std::map<unsigned int, shp<WindkesselModel>>     M_models;
};

}

#endif // BCMANAGER_HPP
