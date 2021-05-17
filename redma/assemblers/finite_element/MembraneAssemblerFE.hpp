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

#ifndef REDMA_MEMBRANEASSEMBLERFE_HPP
#define REDMA_MEMBRANEASSEMBLERFE_HPP

#include <redma/RedMA.hpp>
#include <redma/assemblers/finite_element/NavierStokesAssemblerFE.hpp>
#include <redma/solver/time_marching_algorithms/TimeMarchingAlgorithmFactory.hpp>
#include <redma/geometry/MembraneThicknessComputer.hpp>

namespace RedMA
{

/*! \brief Finite element assembler of the Navier-Stokes problem with the Coupled Momentum method.
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
class MembraneAssemblerFE : public NavierStokesAssemblerFE {

public:
    /*! \brief Constructor taking a datafile, a TreeNode, and a string describing
     * a stabilization as arguments.
     *
     * \param datafile The datafile.
     * \param datafile The TreeNode encoding the physical domain.
     */
    MembraneAssemblerFE(const DataContainer& datafile,
                        shp<TreeNode> treeNode);

    /*! \brief Virtual setup function.
     *
     * Firstly, this calls the parent method in NavierStokesAsseblerFE class.
     * Additionally, it computes and exports the membrane thickness by solving a laplacian
     * problem on the vessel wall, setting the 10% of the tube radius as BC at the extremal rings.
     * Also, it initializes the time-marching scheme employed to extrapolate and update the
     * wall displacement.
     */
    void setup() override;

    /*! \brief Assemble the mass matrix.
     *
     * The mass matrix is defined as
     *
     * \f[
     *    M = M^{NS} + M^{bd}
     * \f]
     *
     * where M_{NS} is the "standard" mass matrix, while M_{bd} is the boundary mass matrix.
     * The former is defined as in \see StokesAssemblerFE::assembleMass(); the latter is defined as
     *
     * \f[
     *    M^{bd}_{ij} = \int_{\Gamma_W} \rho_W H \varphi_j^h \cdot \varphi_i^h
     * \f]
     *
     * where \f$\Gamma_W$\f is the vessel wall, \f$\rho_W$\f is the wall density
     * and \f$H$\f is the (non-constant) wall thickness. \f$\varphi_i^h\f$ are the finite
     * element basis functions of the velocity.
     *
     * \param bcManager A BCManager for the application of the boundary conditions.
     */
    shp<aMatrix> assembleMass(shp<BCManager> bcManager) override;

    /*! \brief Assemble the boundary mass matrix.
     *
     * The boundary mass matrix is defined as
     *
     * \f[
     *    M^{bd}_{ij} = \int_{\Gamma_W} \rho_W H \varphi_j^h \cdot \varphi_i^h
     * \f]
     *
     * where \f$\Gamma_W$\f is the vessel wall, \f$\rho_W$\f is the wall density
     * and \f$H$\f is the (non-constant) wall thickness. \f$\varphi_i^h\f$ are the finite
     * element basis functions of the velocity.
     *
     * \param bcManager A BCManager for the application of the boundary conditions.
     */
    shp<aMatrix> assembleBoundaryMass(shp<BCManager> bcManager, bool verbose = false);

    /*! \brief Assemble the wall boundary mass matrix, i.e. the boundary mass without
     * the contribution of wall density and wall thickness.
     *
     * The wall boundary mass matrix is defined as
     *
     * \f[
     *    M^{wbd}_{ij} = \int_{\Gamma_W} \varphi_j^h \cdot \varphi_i^h
     * \f]
     *
     * where \f$\Gamma_W$\f is the vessel wall. \f$\varphi_i^h\f$ are the finite
     * element basis functions of the velocity.
     *
     * \param bcManager A BCManager for the application of the boundary conditions.
     */
    shp<aMatrix> assembleWallBoundaryMass(shp<BCManager> bcManager, bool verbose = false);

    /*! \brief Assemble the stiffness matrix.
     *
     * The stiffness matrix is defined as
     *
     * \f[
     *    K = K^{NS} + K^{bd}
     * \f]
     *
     * where K_{NS} is the "standard" stiffness matrix, while K_{bd} is the boundary stiffness matrix.
     * The former is defined as in \see StokesAssemblerFE::assembleStiffness(); the latter is defined as
     *
     * \f[
     *    K^{bd}_{ij} = \int_{\Gamma_W} \frac{E}{1+\nu} \frac{\nabla_{\gamma}\varphi_j^h +
     *    \nabla^T_{gamma} \varphi_j^h}{2} : \nabla_{gamma} \varphi_j^h +
     *    \int_{\Gamma_W} \frac{E\nu}{1-\nu^2} (\nabla_{gamma} \cdot \varphi_j^h)(\nabla_{gamma} \cdot \varphi_i^h) +
     *    \int_{\Gamma_W} (k-1)\frac{E}{1+\nu} \frac{\nabla_{\gamma}\varphi_j^h +
     *    \nabla^T_{gamma} \varphi_j^h}{2} n \otimes n : \nabla_{gamma} \varphi_i^h
     * \f]
     *
     * where \f$\Gamma_W$\f is the vessel wall, \f$\nabla_{gamma}$\f denotes the surface gradient,
     * \f$E$\f is the Young modulus, \f$\nu$\f is the Poisson ratio and \f$n$\f is the outward
     * unit normal vector. \f$\varphi_i^h\f$ are the finite element basis functions of the velocity.
     *
     * \param bcManager A BCManager for the application of the boundary conditions.
     */
    shp<aMatrix> assembleStiffness(shp<BCManager> bcManager) override;

    /*! \brief Assemble the boundary stiffness matrix.
     *
     * The boundary stiffness matrix is defined as
     *
     * \f[
     *    K^{bd}_{ij} = \int_{\Gamma_W} \frac{E}{1+\nu} \frac{\nabla_{\gamma}\varphi_j^h +
     *    \nabla^T_{gamma} \varphi_j^h}{2} : \nabla_{gamma} \varphi_j^h +
     *    \int_{\Gamma_W} \frac{E\nu}{1-\nu^2} (\nabla_{gamma} \cdot \varphi_j^h)(\nabla_{gamma} \cdot \varphi_i^h) +
     *    \int_{\Gamma_W} (k-1)\frac{E}{1+\nu} \frac{\nabla_{\gamma}\varphi_j^h +
     *    \nabla^T_{gamma} \varphi_j^h}{2} n \otimes n : \nabla_{gamma} \varphi_i^h
     * \f]
     *
     * where \f$\Gamma_W$\f is the vessel wall, \f$\nabla_{gamma}$\f denotes the surface gradient,
     * \f$E$\f is the Young modulus, \f$\nu$\f is the Poisson ratio and \f$n$\f is the outward
     * unit normal vector. \f$\varphi_i^h\f$ are the finite element basis functions of the velocity.
     *
     * \param bcManager A BCManager for the application of the boundary conditions.
     */
    shp<aMatrix> assembleBoundaryStiffness(shp<BCManager> bcManager, bool verbose = false);

    /*! \brief Getter for the right-hand side term in Newton-Raphson iterations.
     *
     * The right-hand side term for the Newton-Raphson iterations of the Coupled Momentum method
     * is built as follows
     *
     * \f[
     *    rhs(w) = rhs_{NS}(w) + K_{bd} \left(\sum_{j=1}^{M} \alpha_j d_j\right) -
     *             - (\Delta t \beta c_s + k_s) M^{wbd} w +
     *             + c_s M^{wbd} \left(\sum_{j=1}^{M} \alpha_j d_j\right)
     * \f]
     *
     * where \f$\left(\sum_{j=1}^{M} \alpha_j d_j\right)$\f is the combination of the previous
     * wall displacements, according to the chosen time-marching scheme, \f$c_s$\f is the vessel
     * wall elasticity, \f$k_s$\f is the vessel wall visco-elasticity, \f$\beta$\f is the
     * coefficient related to the current time instant in the chosen time-marching scheme.
     * \f$rhs_{NS}$\f is defined as in \see NavierStokesAssemblerFE::getRightHandSide.
     *
     * \param time Current time.
     * \param sol Current solution.
     * \return Shared pointer to aVector of the right-hand side term
     */
    shp<aVector> getRightHandSide(const double& time,
                                  const shp<aVector>& sol) override;

    /*! \brief Getter for Jacobian of the right-hand side term in Newton-Raphson iterations.
     *
     * The right-hand side term Jacobian for the Coupled Momentum method is built as follows:
     *
     * \f[
     *    Jrhs(w) = Jrhs_{NS}(w) - (\Delta t \beta c_s + k_s) M^{wbd}
     * \f]
     *
     * where \f$c_s$\f is the vessel wall elasticity, \f$k_s$\f is the vessel wall visco-elasticity,
     * \f$\beta$\f is the coefficient related to the current time instant in the chosen time-marching
     * scheme. \f$Jrhs_{NS}$\f is defined as in \see NavierStokesAssemblerFE::getJacobianRightHandSide.
     *
     * \param time Current time.
     * \param sol Current solution.
     * \return Shared pointer to aMatrix of the right-hand side term Jacobian.
     */
    shp<aMatrix> getJacobianRightHandSide(const double& time,
                                          const shp<aVector>& sol) override;

    /*! \brief Getter for a vector containing the boundary mass, the boundary stiffness and
    * the wall boundary mass (in order).
    *
    * \return A vector containing the shared pointers to the boundary matrices.
    */
    std::vector<shp<aMatrix>> getBoundaryMatrices() const;

    /*! PostProcess function, to be called at the end of each timestep.
     *
     * The function first calls the parent method \see NavierStokesAssemblerFE::postProcess. Then
     * it also performs the update of the wall displacements, using the same time-marching scheme
     * employed for velocity and pressure, and adds the wall displacement to the exported solution.
     *
     * \param time Current time.
     * \param sol Current solution.
     */
    void postProcess(const double& t, const shp<aVector>& sol) override;

    /*! \brief Setter for the exporter(s).
     *
     * In addition to the Navier-Stokes exporter (\see NavierStokesAssemblerFE::setExporter) which sets
     * exporters for velocity, pressure and wall shear stress (WSS), here also an exporter for the
     * wall displacement is set.
     */
    void setExporter() override;

    /*! \brief Setter method for the wall displacement exporter
     *
     * @param displacement Shared pointer to the EpetraVector storing the wall displacements
     */
    void setDisplacementExporter(shp<VECTOREPETRA> displacement);

    /*! \brief Getter method for the elasticity of the vessel wall
     *
     * @return Elasticity of the vessel wall
     */
    inline double getWallElasticity() const {return M_wall_elasticity;}

    /*! \brief Getter method for the visco-elasticity of the vessel wall
     *
     * @return Visco-elasticity of the vessel wall
     */
    inline double getWallViscoelasticity() const {return M_wall_viscoelasticity;}

    /*! \brief Getter method for the indicator vector of the vessel wall and of the I/O rings
     *
     * @return Shared pointer to an EpetraVector with value of 1 at the vessel wall and I/O rings DOFs
     * and of 0 at all the others
     */
    inline shp<VECTOREPETRA> getBoundaryIndicator() const {return M_boundaryIndicator;}

protected:

    /*! \brief Method that computes the (modified) Lamé constants for the Coupled Momentum method.
     *
     * The (modified) Lamé constants for the Coupled Momentum method are:
     *
     * \f[
     *    \lambda_1 = \frac{E\nu}{1-\nu^2};   \lambda_2 = \frac{E}{2(1+\nu)}
     * \f]
     *
     * where \f$E$\f is the Young modulus and \f$\nu$\f is the Poisson ratio.
     *
     */
    void computeLameConstants();

    /*! \brief Method to compute the membrane thickness.
     *
     * The membrane thickness \f$H$\f is computed by solving a homogeneous Laplacian problem on the
     * vessel wall. The 10% of the vessel diameter is imposed as BC at the extremal I/O rings.
     *
     */
    void computeThickness();

    /*! \brief Method to export the membrane thickness for visualization purposes.
     *
     */
    void exportThickness();

    shp<aTimeMarchingAlgorithm>                     M_TMA_Displacements;

    double                                          M_lameI;
    double                                          M_lameII;
    double                                          M_membrane_density;
    double                                          M_transverse_shear_coeff;
    double                                          M_wall_elasticity;
    double                                          M_wall_viscoelasticity;
    unsigned int                                    M_wallFlag;

    shp<BlockMatrix>                                M_boundaryStiffness;
    shp<BlockMatrix>                                M_boundaryMass;
    shp<BlockMatrix>                                M_wallBoundaryMass;

    shp<VECTOREPETRA>                               M_displacementExporter;
    shp<VECTOREPETRA>                               M_boundaryIndicator;

};

};


#endif //REDMA_MEMBRANEASSEMBLERFE_HPP
