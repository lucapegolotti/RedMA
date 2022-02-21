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

/*! \brief Reduced Basis assembler of the Navier-Stokes problem with the Coupled Momentum method.
 *
 * The equations are defined, for every \f$(x,t)\in \Omega\f$, as
 *
 * \f{eqnarray*}{
 *        \rho \dot{u} + \rho (u \cdot \nabla)u - \nabla \cdot \sigma(u,p) &= f, \\
 *        \nabla \cdot u &= 0,
 * \f}
 * where \f$u\f$ and \f$p\f$ are velocity and pressure of the fluid,
 * \f[
 *       \sigma(u,p) = 2\mu \varepsilon(u) - pI, \quad \varepsilon(u) = (\nabla u + (\nabla u)^T)/2,
 * \f]
 *
 * \f$\rho\f$ is the density of the fluid, and \f$\mu\f$ its viscosity.
 *
 * Under the assumption of small displacements, the Coupled Momentum method allows to embed the
 * dynamics of the vessel motion in the equations for the blood flow. In particular, it allows
 * to obtain a reduced FSI formulation, where a Navier-Stokes system in a fixed fluid domain
 * is supplemented by a Robin boundary condition that acts as a surrogate of the structure model.
 * We refer to the following work for further information:
 *
 * Colciago, C.M., Deparis, S. and Quarteroni, A., 2014.
 * Comparisons between reduced order models and full 3D models for fluid–structure interaction
 * problems in haemodynamics.
 * Journal of Computational and Applied Mathematics, 265, pp.120-138.
 *
 * In particular, the wall displacement is not directly inserted in the system to be solved;
 * rather, a time-marching scheme analogous to the one employed for velocity and pressure
 * is used to perform extrapolation (before the system resolution) and update (after the
 * system resolution).
 *
 * In some methods, it is required to associate variables and matrices to specific
 * indices.
 *
 * Indices of the components: velocity = 0, pressure = 1.
 * Indices of the matrices: mass = 0, stiffness = 1, divergence = 2.
 * Indices of boundary matrices: mass = 0, stiffness = 1, wall mass = 1
 *
 */
class MembraneAssemblerRB : public NavierStokesAssemblerRB {

public:
    /*! \brief Constructor taking a datafile and a TreeNode as arguments.
     *
     * \param datafile The datafile.
     * \param datafile The TreeNode encoding the physical domain.
     */
    MembraneAssemblerRB(const DataContainer &data,
                        shp<TreeNode> treeNode);

    /*! \brief Perform the setup for the RB method.
     *
     * Callback to \see NavierStokesAssemblerFE::RBsetup. Additionally, the boundary matrices
     * (i.e. boundary mass, boundary stiffness, wall boundary mass) employed by the Coupled Momentum
     * method are assembled and projected onto the dimensionality reduced subspace. Finally, the
     * time-marching algorithm responsible for the extrapolation and update of the wall displacements
     * is initialized.
     */
    void RBsetup() override;

    /*! \brief Getter for the right-hand side term.
     *
     * The right-hand side term is assembled as in \see MembraneAssemblerFE::getRightHandSide, except
     * from the fact that all quantities are now projected onto a suitable dimensionality-reduced subspace.
     *
     * \param time Current time.
     * \param sol Current solution.
     * \return Shared pointer to aVector of the right-hand side term
     */
    shp<aVector> getRightHandSide(const double &time,
                                  const shp<aVector> &sol) override;

    /*! \brief VGetter for Jacobian of the right-hand side.
     *
     * Callback to \see NavierStokesAssemblerRB::getJacobianRightHandSide. (Attention: the Jacobian is not
     * consistent as we neglect the convective term). Additionally, suitable terms coming from the Coupled
     * Momentum method are added, as in \see MembraneAssemblerFE::getJacobianRightHandSide, with the
     * difference that here all quantities are projected onto a dimensionality-reduced subspace.
     *
     * \param time Current time.
     * \param sol Current solution.
     * \return Shared pointer to aMatrix of the right-hand side term Jacobian.
     */
    shp<aMatrix> getJacobianRightHandSide(const double& time,
                                          const shp<aVector>& sol,
                                          const double& diagCoeff = 0) override;

    /*! PostProcess function, to be called at the end of each timestep.
     *
     * The function first calls the parent method \see StokesAssemblerRB::postProcess. Then
     * it also performs the update of the wall displacements, using the same time-marching scheme
     * employed for velocity and pressure, and adds the wall displacement to the exported solution.
     *
     * \param time Current time.
     * \param sol Current solution.
     */
    void postProcess(const double& t, const shp<aVector>& sol) override;

protected:

    /*! \brief Getter method for the internal FE Assembler, casted to a MembraneAssemblerFE
     *
     * @return Shared pointer to a MembraneAssemblerFE object, representing the inner FE assembler
     */
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
