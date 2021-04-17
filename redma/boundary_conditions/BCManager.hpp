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

#include <lifev/core/fem/BCHandler.hpp>

namespace RedMA
{

/*! \brief Class to handle the boundary conditions of an assembler.
 */
class BCManager
{
    typedef aTimeMarchingAlgorithm                  TimeMarchingAlgorithm;
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
     */
    void applyDirichletBCs(const double& time,
                           BlockVector& input,
                           shp<FESPACE> fespace,
                           const unsigned int& index) const;

   /*! \brief Apply homogeneous Dirichlet boundary conditions to a block vector.
    *
    * \param input The input block vector.
    * \param fespace The finite element space.
    * \param index Index of the block vector to which the boundary conditions
    * must be applied.
    */
    void apply0DirichletBCs(BlockVector& input,
                            shp<FESPACE> fespace,
                            const unsigned int& index) const;

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
                               const double& diagCoefficient) const;

    /*! \brief Get Neumann condition.
     *
     * The Neumann condition is constant on the corresponding flag.
     *
     * \param time The current time.
     * \param flag The flag of the outlet.
     * \param rate Flow rate at the outlet.
     */
    double getNeumannBc(const double& time,
                        const double& flag,
                        const double& rate);

    /*! \brief Get Jacobian of the Neumann condition.
     *
     * \param time The current time.
     * \param flag The flag of the outlet.
     * \param rate Flow rate at the outlet.
     */
    double getNeumannJacobian(const double& time,
                              const double& flag,
                              const double& rate);

    /// Post process function to call at the end of a timestep.
    void postProcess();
    shp<VECTOREPETRA> computeBoundaryIndicator(shp<FESPACE> fespace) const;

//inline bool useStrongDirichlet() const {return M_strongDirichlet;}


private:
    static double poiseuille(const double& t,
                             const double& x,
                             const double& y,
                             const double& z,
                             const unsigned int& i,
                             const GeometricFace& face,
                             const Law inflow,
                             const double& coefficient);

    static double fZero(const double& t,
                        const double& x,
                        const double& y,
                        const double& z,
                        const unsigned int& i);

    static double constantFunction(const double& t,
                                   const double& x,
                                   const double& y,
                                   const double& z,
                                   const unsigned int& i,
                                   const double& K);

    static double fZero2(double t);

    shp<LifeV::BCHandler> createBCHandler0Dirichlet() const;


    void addInletBC(shp<LifeV::BCHandler> bcs,
                    const Law& law,
                    GeometricFace inlet) const;

    static double fOne(const double& t, const double& x, const double& y,
                       const double& z, const unsigned int& i);

    void parseNeumannData();

    shp<TreeNode>                                    M_treeNode;
    DataContainer                                    M_data;
    std::map<unsigned int, Law>                      M_inflows;
    std::map<unsigned int, Law>                      M_inflowsDt;

    std::vector<unsigned int>                        M_inletFlags;
    unsigned int                                     M_wallFlag;
    unsigned int                                     M_inletRing;
    unsigned int                                     M_outletRing;

    // key is the outlet index (more than one for bifurcations)
    //std::map<unsigned int,shp<WindkesselModel>>      M_models;

    double                                           M_coefficientInflow;
};

}

#endif // BCMANAGER_HPP
