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
#include <redma/boundary_conditions/aBCModel.hpp>
#include <redma/boundary_conditions/Windkessel/WindkesselModel.hpp>
#include <redma/boundary_conditions/Coronary/CoronaryModel.hpp>

#include <lifev/core/fem/BCHandler.hpp>
#include <lifev/core/array/VectorSmall.hpp>
#include <lifev/core/array/MatrixSmall.hpp>

namespace RedMA
{

/*! \brief Class to handle the boundary conditions of an assembler.
 */
class BCManager
{
    typedef LifeV::VectorSmall<3>     Vector3D;
    typedef LifeV::MatrixSmall<3,3>   Matrix3D;

    typedef std::function<double(double)>           Law;

public:
    /*! \brief Constructor accepting a datafile and a TreeNode.
     *
     * \param datafile The DataContainer.
     * \param treeNode Shared pointer to a TreeNode.
     */
    BCManager(const DataContainer& datafile,
              shp<TreeNode> treeNode);

    /*! \brief Apply Dirichlet boundary conditions to a block vector.
     *
     * \param time The current time.
     * \param input The input block vector.
     * \param fespace The finite element space.
     * \param index Index of the block vector to which the boundary conditions
     * must be applied.
     * \param ringOnly True if the BCs have to be imposed only at the rings
     */
    void applyDirichletBCs(const double& time, BlockVector& input,
                           shp<FESPACE> fespace, const unsigned int& index,
                           const bool& ringOnly = false);

   /*! \brief Apply homogeneous Dirichlet boundary conditions to a block vector.
    *
    * \param input The input block vector.
    * \param fespace The finite element space.
    * \param index Index of the block vector to which the boundary conditions
    * must be applied.
    */
    void apply0DirichletBCs(BlockVector& input,
                            shp<FESPACE> fespace,
                            const unsigned int& index,
                            const bool& ringOnly = false);

    /*! \brief Apply Dirichlet boundary conditions to a block vector.
     *
     * \param input The input block vector.
     * \param fespace The finite element space.
     * \param index Index of the row in the matrix to which the boundary conditions
     * must be applied.
     * \param diagCoefficient Coefficient to put in the diagonal of the matrix.
     */
    void apply0DirichletMatrix(BlockMatrix& input,
                               shp<FESPACE> fespace,
                               const unsigned int& index,
                               const double& diagCoefficient,
                               const bool& ringOnly = false);

    void applyInletDirichletBCs(shp<LifeV::BCHandler> bcs, const Law& law, GeometricFace inlet,
                                const bool& zeroFlag = false) const;

    void applyInletNeumannBCs(shp<LifeV::BCHandler> bcs, const Law& law, GeometricFace inlet,
                              const bool& zeroFlag = false) const;

    void applyOutletNeumannBCs(shp<LifeV::BCHandler> bcs,
                               const bool& zeroFlag = true) const;

    double getOutletNeumannBC(const double& time, const double& flag, const double& rate);

    // actually derivative wrt to flowrate
    double getOutletNeumannJacobian(const double& time, const double& flag, const double& rate);

    void postProcess();

    // inline bool useStrongDirichlet() const {return M_strongDirichlet;}

    inline std::string getInletBCType() const {return M_inletBCType;}

    inline std::map<unsigned int, Law> getInletBCs() const {return M_inletBCs;}

    std::vector<unsigned int> getWallFlags(const bool& withRings = true) const;

    shp<VECTOREPETRA> computeBoundaryIndicator(shp<FESPACE> fespace, const std::vector<unsigned int> flags) const;

    static double fZero(const double& t, const double& x, const double& y,
                        const double& z, const unsigned int& i);

    static double fOne(const double& t, const double& x, const double& y,
                       const double& z, const unsigned int& i);

    static double constantFunction(const double& t, const double& x, const double& y,
                                   const double& z, const unsigned int& i,
                                   const double& K);

    shp<VECTOREPETRA> computeRingsIndicator(shp<FESPACE> fespace,
                                            const unsigned int flag = 999,
                                            const bool onlyExtremal = true) const;

private:
    static double poiseuilleInflow(const double& t, const double& x, const double& y,
                                   const double& z, const unsigned int& i,
                                   const GeometricFace& face,
                                   const Law inflow,
                                   const double& coefficient);

    static double neumannLaw(const double& t, const double& x, const double& y,
                             const double& z, const unsigned int& i,
                             Law inflowLaw);

    std::map<unsigned int, Matrix3D> computeRotationMatrices() const;

    void computeGlobalRotationMatrix(shp<FESPACE> fespace);

    void shiftToNormalTangentialCoordSystem(shp<MATRIXEPETRA> mat, shp<VECTOREPETRA> vec,
                                            shp<FESPACE> fespace);

    void shiftToCartesianCoordSystem(shp<MATRIXEPETRA> mat, shp<VECTOREPETRA> vec,
                                     shp<FESPACE> fespace);

    void checkInletLaw();

    shp<LifeV::BCHandler> createBCHandler0Dirichlet(const bool& ringOnly = false) const;

    shp<LifeV::BCHandler> createBCHandler0DirichletRing() const;

    void addInletBC(shp<LifeV::BCHandler> bcs,
                    const Law& law,
                    GeometricFace inlet, 
                    const bool& ringOnly = false,
                    const bool& zeroFlag = false) const;

    void parseOutletBCData();

    shp<TreeNode>                                    M_treeNode;
    DataContainer                                    M_data;

    std::string                                      M_inletBCType;
    std::string                                      M_ringConstraint;
    std::map<unsigned int, Law>                      M_inletBCs;
    std::map<unsigned int, Law>                      M_outletBCs;
    // bool                                             M_strongDirichlet;

    std::vector<unsigned int>                        M_inletFlags;
    std::vector<unsigned int>                        M_outletFlags;
    std::vector<unsigned int>                        M_trueOutletFlags;

    unsigned int                                     M_wallFlag;
    std::vector<unsigned int>                        M_inletRingFlags;
    std::vector<unsigned int>                        M_outletRingFlags;
    std::vector<unsigned int>                        M_trueOutletRingFlags;

    shp<MATRIXEPETRA>                                M_globalRotationMatrix;

    std::map<unsigned int, double>                   M_coefficientsInflow;

    std::map<unsigned int, shp<aBCModel>>            M_models;
};

}

#endif // BCMANAGER_HPP
