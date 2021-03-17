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

#ifndef NAVIERSTOKESASSEMBLERRB_HPP
#define NAVIERSTOKESASSEMBLERRB_HPP

#include <redma/assemblers/abstract/aAssemblerFE.hpp>
#include <redma/assemblers/reduced_basis/StokesAssemblerRB.hpp>
#include <redma/assemblers/abstract/aAssemblerRB.hpp>
#include <redma/reduced_basis/RBBases.hpp>

namespace RedMA
{

/*! \brief Reduced basis assembler of the Navier-Stokes problem.
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
class NavierStokesAssemblerRB : public StokesAssemblerRB
{
public:
    /*! \brief Constructor taking a datafile and a TreeNode as arguments.
     *
     * \param datafile The datafile.
     * \param datafile The TreeNode encoding the physical domain.
     */
    NavierStokesAssemblerRB(const DataContainer& data,
                            shp<TreeNode> treeNode);

    /*! \brief Add Jacobian of convective matrix to an input matrix.
     *
     * \param sol The current solution.
     * \param mat The matrix to be modified.
     */
    void addConvectiveTermJacobian(shp<aVector> sol,
                                   shp<aMatrix> mat);

    /*! \brief Virtual getter for right-hand side.
     *
     * \param time Current time.
     * \param sol Current solution.
     * \return Shared pointer to aVector of the right-hand side
     */
    shp<aVector> getRightHandSide(const double& time,
                                  const shp<aVector>& sol) override;


    /*! \brief Virtual getter for Jacobian of the right-hand side.
     *
     * Attention! The Jacobian is not consistent as we neglect the convective
     * term.
     *
     * \param time Current time.
     * \param sol Current solution.
     * \return Shared pointer to aMatrix of the right-hand side Jacobian.
     */
    shp<aMatrix> getJacobianRightHandSide(const double& time,
                                          const shp<aVector>& sol) override;

    /*! \brief Perform the setup for the RB method.
     *
     * Callback to RBsetup in StokesAssemblerRB. Additionally, if specified in
     * the datafile, we assemble vectors of the form
     *
     *  \f[
     *    \{c(u)_i\}_{jk} = \int_{\Omega} \rho [(\xi_j^h \cdot \nabla)\xi_k^h] \varphi_i^h,
     *  \f]
     *
     * where \f$\rho\f$ is the density of the fluid, \f$\xi_j^h\f$ are the
     * reduced basis functions of the velocity, and \f$\\varphi_i^h\f$ are the
     * finite element basis functions of the velocity.
     */
    virtual void RBsetup() override;

protected:
    std::vector<std::vector<shp<BlockVector>>>       M_nonLinearTermsDecomposition;
    shp<BlockVector>                                 M_nonLinearTerm;
    bool                                             M_exactJacobian;
};

}

#endif // NAVIERSTOKESASSEMBLERRB_HPP
