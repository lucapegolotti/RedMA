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

#ifndef STOKESASSEMBLERRB_HPP
#define STOKESASSEMBLERRB_HPP

#include <redma/assemblers/abstract/aAssemblerFE.hpp>
#include <redma/assemblers/finite_element/StokesAssemblerFE.hpp>
#include <redma/assemblers/abstract/aAssemblerRB.hpp>
#include <redma/reduced_basis/RBBases.hpp>

namespace RedMA
{

/*! \brief Reduced basis assembler of the Stokes problem.
 *
 * The equations are defined, for every \f$(x,t)\in \Omega\f$, as
 *
 * \f{eqnarray*}{
 *        \rho \dot{u}-\nabla \cdot \sigma(u,p) &= f, \\
 *        \nabla \cdot u &= 0,
 * \f}
 * where \f$u\f$ and \f$p\f$ are velocity and pressure of the fluid,
 * \f[
        \sigma(u,p) = 2\mu \varepsilon(u) - pI, \quad \varepsilon(u) = (\nabla u + (\nabla u)^T)/2,
 * \f]
 *
 * \f$\rho\f$ is the density of the fluid, and \f$\mu\f$ its viscosity.
 *
 * In some methods, it is required to associate variables and matrices to specific
 * indices.
 *
 * Indices of the components: velocity = 0, pressure = 1.
 * Indices of the matrices: mass = 0, stiffness = 1, divergence = 2.
 */
class StokesAssemblerRB : public aAssemblerRB
{
public:
    /*! \brief Constructor taking a datafile and a TreeNode as arguments.
     *
     * \param datafile The datafile.
     * \param datafile The TreeNode encoding the physical domain.
     */
    StokesAssemblerRB(const DataContainer& data,
                      shp<TreeNode> treeNode);

    /*! \brief Virtual setup function.
     *
     * In this method, we call the setup method of the internal StokesAssemblerFE.
     */
    virtual void setup() override;

    /*! Virtual export solution.
     *
     * \param time Current time.
     * \param sol Current solution.
     */
    virtual void exportSolution(const double& t,
                                const shp<aVector>& sol) override;

    /*! Virtual postProcess functions (to be called at the end of the timestep).
     *
     * \param time Current time.
     * \param sol Current solution.
     */
    virtual void postProcess(const double& t,
                             const shp<aVector>& sol) override;

    /*! \brief Virtual getter for mass matrix.
     *
     * \param time Current time.
     * \param sol Current solution.
     * \return Shared pointer to aMatrix of the mass matrix.
     */
    virtual shp<aMatrix> getMass(const double& time,
                                 const shp<aVector>& sol) override;

    /*! \brief Virtual getter for mass matrix Jacobian.
     *
     * \param time Current time.
     * \param sol Current solution.
     * \return Shared pointer to aMatrix of the mass matrix Jacobian.
     */
    virtual shp<aMatrix> getMassJacobian(const double& time,
                                         const shp<aVector>& sol) override;

    /*! \brief Virtual getter for right-hand side.
     *
     * \param time Current time.
     * \param sol Current solution.
     * \return Shared pointer to aVector of the right-hand side
     */
    virtual shp<aVector> getRightHandSide(const double& time,
                                          const shp<aVector>& sol) override;

    /*! \brief Virtual getter for Jacobian of the right-hand side.
     *
     * \param time Current time.
     * \param sol Current solution.
     * \return Shared pointer to aMatrix of the right-hand side Jacobian.
     */
    virtual shp<aMatrix> getJacobianRightHandSide(const double& time,
                                                  const shp<aVector>& sol) override;

    /*! \brief Virtual getter for the zero vector.
     *
     * \return Shared pointer to aVector of zeros.
     */
    virtual shp<aVector> getZeroVector() const override;

    /*! \brief Virtual getter for the (finite element) lifting.
     *
     * Callback ti getLifting in the internal StokesAssemblerFE.
     *
     * \param time Current time.
     * \return Shared pointer to aVector of the lifting.
     */
    virtual shp<aVector> getLifting(const double& time) const override;

    /*! \brief Getter for the finite element space corresponding to the Dirichlet bcs,
     * i.e., the velocity finite element space.
     *
     * \return Shared pointer to desired finite element space.
     */
    virtual inline shp<FESPACE> getFESpaceBCs() const override
    {
        return M_feStokesAssembler->getFESpaceBCs();
    }

    /*! \brief Getter for the component associated with the Dirichlet bcs.
     *
     * \return Index of the component associated with the Dirichlet bcs (namely index of velocity, 0).
     */
    virtual inline unsigned int getComponentBCs() const override {return 0;}

    /*! \brief Getter for the ETA finite element space associated with the coupling (velocity).
     *
     * \return Shared pointer to desired ETA finite element space.
     */
    virtual inline shp<ETFESPACE3> getETFESpaceCoupling() const override
    {
        return M_feStokesAssembler->getETFESpaceCoupling();
    }

    /*! \brief Virtual method to apply Dirichlet bcs to a matrix.
     *
     * Empty method.
     *
     * \param matrix The matrix to which the bcs must be applied.
     * \param diagCoeff Coefficient to put in the diagonal of the matrix.
     */
    void applyDirichletBCsMatrix(shp<aMatrix> matrix,
                                 double diagCoeff) const override;

    /*! \brief Virtual method to apply homogeneous Dirichlet bcs to a vector.
     *
     * Empty method.
     *
     * \param vector The vector to which the bcs must be applied.
     */
    void apply0DirichletBCs(shp<aVector> vector) const override;

    /*! \brief Virtual method to apply Dirichlet bcs to a vector.
     *
     * Empty method.
     *
     * \param time Current time.
     * \param vector The vector to which the bcs must be applied.
     */
    void applyDirichletBCs(const double& time,
                           shp<aVector> vector) const override;

    /*! \brief Getter for the finite element space corresponding to a specific component.
     *
     * \param index Index of the desired component.
     * \return Shared pointer to desired finite element space.
     */
    virtual inline shp<FESPACE> getFEspace(unsigned int index) const override
    {
        return M_feStokesAssembler->getFEspace(index);
    }

    /*! \brief Getter for a vector containing the mass, stiffness and divergence matrix
     *  (in order).
     *
     * \return A vector containing the shared pointers to the matrices.
     */
    virtual std::vector<shp<aMatrix>> getMatrices() const override;

    /*! \brief Assemble matrix corresponding to a specific index.
     *
     * Empty method.
     *
     * \param index Index of matrix. Mass = 0, stiffness = 1, divergence = 2.
     * \return Shared pointer to the matrix.
     */
    virtual shp<aMatrix> assembleMatrix(const unsigned int& index) override;

    /*! \brief Set extrapolated solution.
     *
     * Method not implemented (raises an exception).
     *
     * \param exSol The extrapolated solution.
     */
    void setExtrapolatedSolution(const shp<aVector>& exSol) override
    {
        throw new Exception("setExtrapolatedSolution method not implemented for RB");
    }

    /*! \brief Applies the piola transformation (or its inverse) to a vector.
     *
     * Callback to applyPiola in StokesAssemblerFE.
     *
     * \param solution Shared pointer to aVector to transform.
     * \param inverse If true, inverse of Piola transformation is applied.
     */
    virtual void applyPiola(shp<aVector> solution,
                            bool inverse) override;

    /*! \brief Perform setup for RB method.
     *
     * In this method: i) we scale the velocity vectors by the Piola transformation,
     * and ii) we project the finite element matrices onto the RB space.
     */
    virtual void RBsetup() override;

    /*! \brief Getter for the reduced bases.
     *
     * \return Shared pointer to the RBBases.
     */
    virtual shp<RBBases> getRBBases() const override;

    /*! \brief Setter for the reduced bases.
     *
     * The method accepts an RBBasesManager, and the correct RBBases are selected
     * according to the mesh of the current block.
     *
     * \param rbManager Shared pointer to a RBBasesManager.
     */
    virtual void setRBBases(shp<RBBasesManager> rbManager) override;

    /*! \brief Retrieve finite element vector from a reduced basis vector.
     *
     * \param rbSolution Shared pointer to reduced basis dofs.
     * \return Shared pointer to finite element vector.
     */
    virtual shp<aVector> convertFunctionRBtoFEM(shp<aVector> rbSolution) const override;

    /* \brief Setter for the default assemblers.
     *
     * \param Shared pointer to the DefaultAssemblersLibrary.
     */
    virtual void setDefaultAssemblers(shp<DefaultAssemblersLibrary> defAssemblers) override;

    /*! \brief Getter for the BCManager of the internal StokesAssemblerFE.
     *
     * \return Shared pointer to the BCManager.
     */
    inline virtual shp<BCManager> getBCManager() const override
    {
        return M_feStokesAssembler->getBCManager();
    }

protected:
    shp<LifeV::Exporter<MESH>>                        M_exporter;
    shp<VECTOREPETRA>                                 M_velocityExporter;
    shp<VECTOREPETRA>                                 M_WSSExporter;
    shp<VECTOREPETRA>                                 M_pressureExporter;
    std::string                                       M_name;
    shp<BlockVector>                                  M_extrapolatedSolution;
    shp<RBBases>                                      M_bases;
    shp<BlockMatrix>                                  M_reducedMass;
    shp<BlockMatrix>                                  M_reducedDivergence;
    shp<BlockMatrix>                                  M_reducedStiffness;
    shp<StokesAssemblerFE>                            M_feStokesAssembler;
};

}

#endif // STOKESASSEMBLERRB_HPP
