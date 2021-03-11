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

#ifndef NAVIERSTOKESASSEMBLERFE_HPP
#define NAVIERSTOKESASSEMBLERFE_HPP

#include <redma/assemblers/finite_element/StokesAssemblerFE.hpp>
#include <redma/assemblers/finite_element/SUPGStabilization.hpp>

namespace RedMA
{

/*! \brief Finite element assembler of the Navier-Stokes problem.
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
 * In some methods, it is required to associate variables and matrices to specific
 * indices.
 *
 * Indices of the components: velocity = 0, pressure = 1.
 * Indices of the matrices: mass = 0, stiffness = 1, divergence = 2.
 *
 */
class NavierStokesAssemblerFE : public StokesAssemblerFE
{
public:
    /*! \brief Constructor taking a datafile, a TreeNode, and a string describing
     * a stabilization  as arguments.
     *
     * \param datafile The datafile.
     * \param datafile The TreeNode encoding the physical domain.
     * \param stabilizationName Name of the stabilization.
     */
    NavierStokesAssemblerFE(const DataContainer& data, shp<TreeNode> treeNode,
                            std::string stabilizationName = "");

    /*! \brief Virtual setup function.
     *
     * This is basically a callback to the setup function in StokesAssemblerFE.
     * Additionally, we initialize the structures for the stabilization if needed.
     */
    virtual void setup() override;

    /*! \brief Add convective matrix to an input matrix.
     *
     * The convective matrix is defined as
     *  \f[
     *    C_{ij}(u) = \int_{\Omega} \rho [(u \cdot \nabla)\varphi_j^h] \varphi_i^h,
     *  \f]
     *
     * where \f$\rho\f$ is the fluid density and \f$\varphi_i^h\f$ are the finite
     * element basis functions of the velocity.
     *
     * \param sol The current solution.
     * \param mat The matrix to be modified.
     */
    void addConvectiveMatrix(shp<aVector> sol, shp<aMatrix> mat);


    /*! \brief Add Jacobian of convective matrix to an input matrix.
     *
     * \param sol The current solution.
     * \param mat The matrix to be modified.
     */
    void addConvectiveTermJacobianRightHandSide(shp<aVector> sol,
                                                shp<aMatrix> mat);

    /*! \brief Virtual getter for mass matrix.
     *
     * See StokesAssemblerFE for the definition of the mass matrix.
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


protected:
    shp<NavierStokesStabilization>                    M_stabilization;
    std::string                                       M_stabilizationName;
};

}

#endif // NAVIERSTOKESASSEMBLERFE_HPP
