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

#ifndef aASSEMBLER_HPP
#define aASSEMBLER_HPP

#include <redma/RedMA.hpp>
#include <redma/array/BlockMatrix.hpp>
#include <redma/boundary_conditions/BCManager.hpp>
#include <redma/solver/time_marching_algorithms/aFunctionProvider.hpp>
#include <redma/geometry/TreeStructure.hpp>

#include <redma/reduced_basis/RBBasesManager.hpp>

#include <redma/problem/DataContainer.hpp>
#include <redma/assemblers/DefaultAssemblersLibrary.hpp>

namespace RedMA
{

/*! \brief Abstract assembler class.
 *
 * This class takes care of the assembly of the structures related to a particular
 * discretized PDE. We recall that we consider PDEs of the form
 *  \f[
 *    M(u) \dot{u} = F(t,u),
 *  \f]
 * where \f$M\f$ is the mass matrix, \f$u\f$ is the solution and \f$F\f$ is the right-hand side.
 */
class aAssembler : public aFunctionProvider
{
    typedef DefaultAssemblersLibrary DefaultAssemblers;
public:
    /// \brief Default (empty) constructor
    aAssembler() {}

    /*! \brief Constructor taking a datafile as argument.
     *
     * \param datafile The datafile.
     */
    aAssembler(const DataContainer& datafile);

    /*! \brief Constructor taking a datafile and a TreeNode as argument.
     *
     * \param datafile The datafile.
     * \param datafile The TreeNode encoding the physical domain.
     */
    aAssembler(const DataContainer& datafile,
               shp<TreeNode> node);

    /// Virtual setup function.
    virtual void setup() = 0;

    /*! Virtual export solution.
     *
     * \param time Current time.
     * \param sol Current solution.
     */
    virtual void exportSolution(const double& time,
                                const shp<aVector>& sol) = 0;

    /*! Virtual postProcess functions (to be called at the end of the timestep).
     *
     * \param t Current time.
     * \param sol Current solution.
     */
    virtual void postProcess(const double& t,
                             const shp<aVector>& sol) = 0;

    /*! \brief Virtual getter for mass matrix.
     *
     * \param time Current time.
     * \param sol Current solution.
     * \return Shared pointer to aMatrix of the mass matrix.
     */
    virtual shp<aMatrix> getMass(const double& time,
                                 const shp<aVector>& sol) = 0;

    /*! \brief Virtual getter for mass matrix jacobian.
     *
     * \param time Current time.
     * \param sol Current solution.
     * \return Shared pointer to aMatrix of the mass matrix jacobian.
     */
    virtual shp<aMatrix> getMassJacobian(const double& time,
                                         const shp<aVector>& sol) = 0;

    /*! \brief Virtual getter for right-hand side.
     *
     * \param time Current time.
     * \param sol Current solution.
     * \return Shared pointer to aVector of the right-hand side
     */
    virtual shp<aVector> getRightHandSide(const double& time,
                                          const shp<aVector>& sol) = 0;

    /*! \brief Virtual getter for Jacobian of the right-hand side.
     *
     * \param time Current time.
     * \param sol Current solution.
     * \return Shared pointer to aMatrix of the right-hand side jacobian.
     */
    virtual shp<aMatrix> getJacobianRightHandSide(const double& time,
                                                  const shp<aVector>& sol) = 0;

    /*! \brief Virtual getter for the lifting.
     *
     * \param time Current time.
     * \return Shared pointer to aVector of the lifting.
     */
    virtual shp<aVector> getLifting(const double& time) const = 0;

    /*! \brief Virtual getter for the zero vector.
     *
     * \return Shared pointer to aVector of zeros.
     */
    virtual shp<aVector> getZeroVector() const = 0;

    /*! Getter for the finite element space corresponding to a specific component.
     *
     * \param index Index of the desired component.
     * \return Shared pointer to desired finite element space.
     */
    virtual inline shp<FESPACE> getFEspace(unsigned int index) const {return nullptr;}

    /*! \brief Getter for the finite element space corresponding to the Dirichlet bcs (e.g., in the
     *         Stokes equations, the finite element s of the velocity).
     *
     * \return Shared pointer to desired finite element space.
     */
    virtual inline shp<FESPACE> getFESpaceBCs() const {return nullptr;}

    /*! \brief Getter for the component associated with the Dirichlet bcs.
     *
     * \return Index of the component associated with the Dirichlet bcs.
     */
    virtual inline unsigned int getComponentBCs() const {return 0;}

    /*! \brief Getter for the internal TreeNode.
     *
     * \return Shared pointer to the TreeNode.
     */
    virtual inline shp<TreeNode> getTreeNode() const {return M_treeNode;}

    /*! Getter for the number of components.
     *
     * \return Number of components in the PDE.
     */
    virtual inline unsigned int getNumComponents() const {return M_nComponents;}

    /*! Getter for the internal MPI Communicator.
     *
     * \return Internal MPI Communicator.
     */
    virtual inline EPETRACOMM getComm() const {return M_comm;}

    /*! \brief Getter for the ETFESpace associated with the coupling.
     *
     * Attention: here we assume that the coupling involves three components.
     * This is not true for scalar PDEs (e.g., Laplacian).
     *
     * \return Shared pointer to the ETFESpace associated with the coupling.
     */
    virtual inline shp<ETFESPACE3> getETFESpaceCoupling() const {return nullptr;}

    /*! \brief Getter for the vector of assembled matrices.
     *
     * If not overloaded, this method returns an empty vector of aMatrix.
     *
     * \return Vector of shared pointers to aMatrix (empty if not overloaded).
     */
    virtual inline std::vector<shp<aMatrix>> getMatrices() const {return std::vector<shp<aMatrix>>();}

    /*! \brief Getter for the internal BCManager.
     *
     * \return Shared pointer to the BCManager.
     */
    virtual inline shp<BCManager> getBCManager() const {return M_bcManager;}

    /*! \brief Virtual setup of the exporter.
     */
    virtual void setExporter() = 0;

    /*! \brief Virtual method to apply Dirichlet bcs to a matrix.
     *
     * \param matrix The matrix to which the bcs must be applied.
     * \param diagCoeff Coefficient to put in the diagonal of the matrix.
     */
    virtual void applyDirichletBCsMatrix(shp<aMatrix> matrix,
                                         double diagCoeff) const = 0;

    /*! \brief Virtual method to apply homogeneous Dirichlet bcs to a vector.
     *
     * \param vector The vector to which the bcs must be applied.
     */
    virtual void apply0DirichletBCs(shp<aVector> vector) const = 0;

    /*! \brief Virtual method to apply Dirichlet bcs to a vector.
     *
     * The boundary condtion is evaluated at the time provided as input.
     *
     * \param time Current time.
     * \param vector The vector to which the bcs must be applied.
     */
    virtual void applyDirichletBCs(const double& time,
                                   shp<aVector> vector) const = 0;

    /*! \brief Applies the piola transformation (or its inverse) to a vector.
     *
     * This method must raise an Exception if the Piola transformation does not
     * make sense for the current problem.
     * \param solution Shared pointer to aVector to transform.
     * \param inverse If true, inverse of Piola transformation is applied.
     */
    virtual void applyPiola(shp<aVector> solution,
                            bool inverse) = 0;

    /*! \brief Assemble matrix corresponding to specific index.
     *
     * \param index Desired index.
     * \param structure BlockMDEIMStructure instance (currently not supported).
     * \return The assembled matrix.
     */
    virtual shp<aMatrix> assembleMatrix(const unsigned int& index) {return shp<aMatrix>();}

    /*! \brief Get nonlinear part of the right-hand side (when applicable).
     *
     * \return Shared pointer to aVector of the nonlinear term.
     */
    virtual shp<aVector> getNonLinearTerm() {};

    /*! \brief Initializer for the finite element spaces.
     *
     * Not implemented for aAssembler. It must be overloaded by the local finite
     * element assemblers (but not the block one).
     */
    virtual void initializeFEspaces() {};

    /*! \brief Setter for the default assemblers.
     *
     * \param Shared pointer to the DefaultAssemblersLibrary.
     */
    virtual void setDefaultAssemblers(shp<DefaultAssemblers> defAssemblers)
    {
        M_defaultAssemblers = defAssemblers;
    };

    /*! \brief Get the ID of the internal TreeNode.
     *
     * \return Index of the tree node.
     */
    inline unsigned int ID() {return M_treeNode->M_ID;}

    /// Perform the setup necessary to RB method.
    virtual void RBsetup() {}

    /*! \brief Getter for the RBBases.
     *
     * \return Shared pointer to RBBases.
     */
    virtual shp<RBBases> getRBBases() const {return nullptr;}

    /*! \brief Setter for the RBBasesManager.
     *
     * \param rbManager Shared pointer to RBBasesManager.
     */
    virtual void setRBBases(shp<RBBasesManager> rbManager) {}

    /*! \brief Getter for the no-slip BCs at the lateral wall.
     *
     * \return True if no-slip BC at the lateral wall are imposed.
     */
    // virtual inline bool hasNoSlipBCs() const = 0;

protected:
    DataContainer                           M_data;
    shp<TreeNode>                           M_treeNode;
    shp<BCManager>                          M_bcManager;
    unsigned int                            M_nComponents;
    EPETRACOMM                              M_comm;
    std::string                             M_name;
    shp<DefaultAssemblers>                  M_defaultAssemblers;
};

}

#endif // aASSEMBLER_HPP
