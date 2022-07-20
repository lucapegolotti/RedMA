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

#ifndef STOKESASSEMBLERFE_HPP
#define STOKESASSEMBLERFE_HPP

#include <redma/RedMA.hpp>
#include <redma/assemblers/abstract/aAssemblerFE.hpp>
#include <redma/assemblers/coupling/InterfaceAssembler.hpp>
#include <lifev/eta/expression/Integrate.hpp>
#include <redma/coupling_basis_functions/BasisFunctionFunctor.hpp>
#include <redma/coupling_basis_functions/weightFunction.hpp>
#include <redma/reduced_basis/RBBases.hpp>

namespace RedMA
{
        typedef LifeV::VectorSmall<3> Vector3D;
        /*! \brief Finite element assembler of the Stokes problem.
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

class StokesAssemblerFE : public aAssemblerFE
{
public:
    /*! \brief Constructor taking a datafile and a TreeNode as arguments.
     *
     * \param datafile The datafile.
     * \param datafile The TreeNode encoding the physical domain.
     */
    StokesAssemblerFE(const DataContainer& data,
                      shp<TreeNode> treeNode);

    /// Virtual setup function.
    virtual void setup() override;

    /*! \brief Solutions importing method
     *
     * Method to import a solutions (velocity or pressure) from txt file into an EpetraMartrix
     */
    std::map<unsigned int, std::vector<shp<aVector>>> importSolution(const std::string& filename) const override;

    /*! Virtual export solution.
     *
     * \param time Current time.
     * \param sol Current solution.
     */
    virtual void exportSolution(const double& time,
                                const shp<aVector>& sol) override;

    /*! \brief Virtual export solution to .txt
    *
    * \param time Current time.
    * \param sol Current solution.
    */
    virtual void exportSolutionToTxt(const double& time,
                                     const shp<aVector>& sol,
                                     const std::string& filename) override;

    /*! Virtual postProcess functions (to be called at the end of the timestep).
     *
     * \param time Current time.
     * \param sol Current solution.
     */
    virtual void postProcess(const double& time,
                             const shp<aVector>& sol) override;

    /*! \brief Virtual getter for mass matrix.
     *
     * \param time Current time.
     * \param sol Current solution.
     * \return Shared pointer to aMatrix of the mass matrix.
     */
    virtual shp<aMatrix> getMass(const double& time,
                                 const shp<aVector>& sol) override;

    /*! \brief Virtual getter for pressure mass matrix.
     *
     * \param time Current time.
     * \param sol Current solution.
     * \return Shared pointer to aMatrix of the pressure mass matrix.
     */
    virtual shp<aMatrix> getPressureMass(const double& time,
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
                                                  const shp<aVector>& sol,
                                                  const double& diagCoeff = 0) override;

    /*! \brief Virtual getter for the zero vector.
     *
     * \return Shared pointer to aVector of zeros.
     */
    virtual shp<aVector> getZeroVector() const override;

    /*! \brief Virtual getter for the lifting.
     *
     * \param time Current time.
     * \return Shared pointer to aVector of the lifting.
     */
    virtual shp<aVector> getLifting(const double& time) const override;

    /*! \brief Getter for the FE lifting.
     *
     * \param time Current time.
     * \return Shared pointer to aVector of the FE lifting.
     */
    virtual shp<aVector> getFELifting(const double& time) const override;

    /*! \brief Initializer for the finite element spaces.
     *
     * Here we create the finite element spaces for velocity and pressure.
     */
    void initializeFEspaces() override;

    /*! \brief Setter for the exporter(s).
     *
     * We define the exporters for velocity, pressure and wall shear stress.
     */
    void setExporter() override;

    /*! \brief Getter for the finite element space corresponding to the Dirichlet bcs,
     * i.e., the velocity finite element space.
     *
     * \return Shared pointer to the desired finite element space.
     */
    virtual inline shp<FESPACE> getFESpaceBCs() const override
    {
        return this->M_velocityFESpace;
    }

    /*! \brief Getter for the component associated with the Dirichlet bcs.
     *
     * \return Index of the component associated with the Dirichlet bcs (namely index of velocity, 0).
     */
    virtual inline unsigned int getComponentBCs() const override {return 0;}

    /*! \brief Getter for the flag telling if no-slip BCs on the lateral wall are imposed
    *
    * \return true if no-slip BCs are imposed on th alteral wall, false otherwise
    */
    virtual inline bool hasNoSlipBCs() const {return M_addNoSlipBC;}

    /*! \brief Getter for the ETA finite element space associated with the coupling (velocity).
     *
     * \return Shared pointer to desired ETA finite element space.
     */
    virtual inline shp<ETFESPACE3> getETFESpaceCoupling() const override
    {
        return this->M_velocityFESpaceETA;
    }

    /*! \brief Virtual method to apply Dirichlet bcs to a matrix.
     *
     * \param matrix The matrix to which the bcs must be applied.
     * \param diagCoeff Coefficient to put in the diagonal of the matrix.
     */
    virtual void applyDirichletBCsMatrix(shp<aMatrix> matrix,
                                         double diagCoeff) const override;

    /*! \brief Virtual method to apply homogeneous Dirichlet bcs to a vector.
     *
     * \param vector The vector to which the bcs must be applied.
     */
    virtual void apply0DirichletBCs(shp<aVector> vector) const override;

    /*! \brief Virtual method to apply Dirichlet bcs to a vector.
     *
     * The boundary condtion is evaluated at the time provided as input.
     *
     * \param time Current time.
     * \param vector The vector to which the bcs must be applied.
     */
    virtual void applyDirichletBCs(const double& time,
                                   shp<aVector> vector) const override;

    /*! \brief Getter for the finite element space corresponding to a specific component.
     *
     * \param index Index of the desired component.
     * \return Shared pointer to desired finite element space.
     */
    virtual shp<FESPACE> getFEspace(unsigned int index) const override;

    /*! \brief Getter for the ETA finite element space associated with velocity.
     *
     * \return Shared pointer to the desired ETA finite element space.
     */
    shp<ETFESPACE3> getVelocityETFEspace() const {return M_velocityFESpaceETA;}

    /*! \brief Getter for the ETA finite element space associated with pressure.
    *
    * \return Shared pointer to the desired ETA finite element space.
    */
    shp<ETFESPACE1> getPressureETFEspace() const {return M_pressureFESpaceETA;}

    /*! \brief Getter for a vector containing the mass, stiffness and divergence matrix
     *  (in order).
     *
     * \return A vector containing the shared pointers to the matrices.
     */
    virtual std::vector<shp<aMatrix>> getMatrices() const override;

    /*! \brief Assemble matrix corresponding to a specific index.
     *
     * \param index Index of matrix. Mass = 0, stiffness = 1, divergence = 2.
     * \return Shared pointer to the matrix.
     */
    virtual shp<aMatrix> assembleMatrix(const unsigned int& index) override;

    /*! \brief Getter for the norm matrices.
     *
     * The matrices are assembled within this function. If fieldIndex = 0 and bcs = true,
     * the Dirichlet boundary conditions are applied to the velocity matrix.
     *
     * \param fieldIndex Index of norm matrix (velocity = 0, pressure = 1).
     * \param bcs If true, Dirichlet boundary conditions are applied (only to velocity).
     * \return Shared pointer to the norm matrix.
     */
    virtual shp<aMatrix> getNorm(const unsigned int& fieldIndex,
                                 bool bcs = true) override;

    /*! \brief Getter for the constraint matrix.
     *
     * In this case, we return the divergence matrix. This method is called,
     * for example, when we need to compute the supremizers for the velocity.
     *
     * \return Shared pointer to the matrix.
     */
    virtual shp<aMatrix> getConstraintMatrix() override;

    /*! \brief Set extrapolated solution.
     *
     * \param exSol The extrapolated solution.
     */
    void setExtrapolatedSolution(const shp<aVector>& exSol) override;

    /*! \brief Applies the piola transformation (or its inverse) to a vector.
     *
     * \param solution Shared pointer to aVector to transform.
     * \param inverse If true, inverse of Piola transformation is applied.
     */
    virtual void applyPiola(shp<aVector> solution,
                            bool inverse) override;

    /*! \brief Add Neumann BCs to a vector.
     *
     * The conditions can depend on time and the solution.
     *
     * \param time The current time.
     * \param sol The current solution.
     * \param rhs The vector to which the boundary conditions must be applied.
     */
    void addNeumannBCs(double time, shp<aVector> sol, shp<aVector> rhs);

    /*! \brief Setter for the velocity finite element order (e.g., "P2")
     *
     * \param The velocity finite element order.
     */
    void setVelocityOrder(std::string velocityOrder) {M_velocityOrder = velocityOrder;}

    /*! \brief Setter for the pressure finite element order (e.g., "P1")
     *
     * \param The pressure finite element order.
     */
    void setPressureOrder(std::string pressureOrder) {M_pressureOrder = pressureOrder;}

    /*! \brief Assemble the stiffness matrix.
     *
     * Depending on the value of fluid/use_strain in the datafile, the
     * stiffness matrix is defined as (use_strain = true)
     *  \f[
     *    K_{ij} = \int_{\Omega} \mu / 2 (\nabla \varphi_j^h + (\nabla \varphi_j^h)^T) : (\nabla \varphi_i^h + (\nabla \varphi_i^h)^T),
     *  \f]
     *
     * or (use_strain = false)
     *
     *  \f[
     *    K_{ij} = \int_{\Omega} \mu \nabla \varphi_j^h : \nabla \varphi_i^h,
     *  \f]
     *
     * where \f$\mu\f$ is the fluid viscosity and \f$\varphi_i^h\f$ are the finite element basis functions of the velocity.
     *
     * \param bcManager A BCManager for the application of the boundary conditions.
     */
    virtual shp<aMatrix> assembleStiffness(shp<BCManager> bcManager);

    /*! \brief Assemble the mass matrix.
     *
     * The mass matrix is defined as
     *
     * \f[
     *    M_{ij} = \int_{\Omega} \rho \varphi_j^h \cdot \varphi_i^h
     * \f]
     *
     * where \f$\rho\f$ is the fluid density and \f$\varphi_i^h\f$ are the finite
     * element basis functions of the velocity.
     *
     * \param bcManager A BCManager for the application of the boundary conditions.
     */
    virtual shp<aMatrix> assembleMass(shp<BCManager> bcManager);

    /*! \brief Assemble the mass matrix for pressure.
     *
     * The mass matrix for pressure is defined as
     *
     * \f[
     *    M_{ij} = \int_{\Omega} \psi_j^h \psi_i^h
     * \f]
     *
     * where \f$\psi_i^h\f$ are the finite element basis functions of the pressure.
     *
     * \param bcManager A BCManager for the application of the boundary conditions.
     */
    virtual shp<aMatrix> assemblePressureMass(shp<BCManager> bcManager);

    /*! \brief Assemble the divergence matrix.
     *
     * The divergence matrix is defined as
     *
     * \f[
     *    D_{ij} = -\int_{\Omega} \psi_j^h \nabla \cdot \varphi_i^h
     * \f]
     *
     * where \f$\psi_i^h\f$ and \f$\varphi_i^h\f$ are the finite element functions of pressure and velocity, respectively.
     * \param bcManager A BCManager for the application of the boundary conditions.
     */
    shp<aMatrix> assembleDivergence(shp<BCManager> bcManager);

    /*! \brief Compute flow rates at inflow and outflows.
     *
     * \param sol Current solution.
     * \param verbose If true, flow rates are printed to terminal.
     * \return A map with key = face flag and value = flowrate.
     */
    virtual std::map<unsigned int, double> computeFlowRates(shp<aVector> sol, bool verbose = false) override;

    /*! \brief Assemble vectors to compute the flow rate.
     *
     * See assembleFlowRateVector for the definition of these vectors.
     */
    void assembleFlowRateVectors();

    /*! \brief Compute the jacobians of the vectors computed in computeFlowRates.
     *
     * Attention: this function is expensive.
     */
    void assembleFlowRateJacobians();

    /*! \brief Compute additional outlet matrices.
     *
     */
    void assembleAdditionalOutletMatrices();

    /*! \brief Assemble vector to compute the flow rate given a face.
     *
     * These are defined as
     *
     * \f[
     *     a_i = \int_{\Gamma} \varphi_i^h \cdot n
     * \f]
     *
     * where \f$\varphi_i^h\f$ is the finite element functions of velocity, and \f$n\f$ is the normal to the face \f$\Gamma\f$.
     *
     * \param face The face.
     * \return Shared pointer to the vector.
     */
    shp<VECTOREPETRA> assembleFlowRateVector(const GeometricFace& face);

    /*! \brief Assemble jacobian of the vector computed in assembleFlowRateVector.
     *
     * \param face The face.
     * \return Shared pointer to the vector.
     */
    shp<MATRIXEPETRA> assembleFlowRateJacobian(const GeometricFace& face);

    /*! \brief Assemble the additional matrix appearing in outlet if non-standard BCs are imposed.
     *
     * The additional outlet matrix is defined as
     *
     * \f[
     *    O_{ij} = -\int_{\Gamma} \mu \psi_j^h \cdot \nabla \varphi_i^h n
     * \f]
     *
     * where \f$\varphi_i^h\f$ is the finite element functions of velocity, and \f$n\f$ is the normal to the face \f$\Gamma\f$.
     *
     * \param face The face.
     * \return Shared pointer to the matrix.
     */
    shp<MATRIXEPETRA> assembleAdditionalOutletMatrix(const GeometricFace& face);

    /*! \brief Add backflow stabilization.
     *
     * Currently not implemented.
     *
     * \param input Vector that must be modified.
     * \param sol Current solution
     * \param faceFlag Flag of the desired face.
     */
    void addBackFlowStabilization(shp<aVector>& input,
                                  shp<aVector> sol,
                                  const unsigned int& faceFlag);

    /*! \brief Exports norm of velocity and pressure to file.
     *
     * This function exports only if exporter/exportnorms in the datafile is true.
     *
     * \param time Current time.
     * \param velocity Current velocity.
     * \param pressure Current pressure.
     */
    void exportNorms(double time, shp<VECTOREPETRA> velocity, shp<VECTOREPETRA> pressure);

    /*! \brief Compute the wall shear stress given the velocity.
     *
     * \param velocity The current velocity.
     * \param WSS Output parameter: the wall shear stress.
     * \param comm The MPI communicator.
     */
    void computeWallShearStress(shp<VECTOREPETRA> velocity,
                                shp<VECTOREPETRA> WSS,
                                EPETRACOMM comm);

    /*! \brief Integrate the wall shear stress on area of interest.
     *
     * \param velocity The current velocity.
     * \param WSS Output parameter: the wall shear stress.
     * \param comm The MPI communicator.
     */
    void integrateWallShearStress(shp<VECTOREPETRA> velocity, shp<VECTOREPETRA> WSS, EPETRACOMM comm, double intWSS);

    /*! \brief weight function for integrateWSS
     *
     * @param x x
     * @param y y
     * @param z z
     * @param center center of integration
     * @param radius radius of integration
     * @return
     */
    static double weightFunction(const double& t,
                          const double& x,
                          const double& y,
                          const double& z,
                          const LifeV::ID& i);

    void exportIntWSS(double intWSS);

    /*! \brief Initialize the finite element space of the velocity.
     *
     * \param comm The MPI communicator.
     */
    void initializeVelocityFESpace(EPETRACOMM comm);

    /*! \brief Initialize the finite element space of the pressure.
     *
     * \param comm The MPI communicator.
     */
    void initializePressureFESpace(EPETRACOMM comm);

    void initializeDisplacementFESpace(EPETRACOMM comm);

    /*! \brief Get the forcing term \f$f\f$
     *
     * \param time The current time.
     * \return The forcing term.
     */
    shp<aVector> getForcingTerm(const double& time) const;

    /*! \brief Getter for the density.
     *
     * \return The density \f$\rho\f$.
     */
    inline double getDensity() {return M_density;}

    /*! \brief Getter for the viscosity.
     *
     * \return The viscosity \f$\mu\f$.
     */
    inline double getViscosity() {return M_viscosity;}

protected:
    shp<BlockVector> buildZeroVector() const;

    shp<LifeV::Exporter<MESH>>                        M_exporter;
    shp<VECTOREPETRA>                                 M_velocityExporter;
    shp<VECTOREPETRA>                                 M_WSSExporter;
    shp<VECTOREPETRA>                                 M_displacementExporter;
    double                                            M_intWSSExporter;
    shp<VECTOREPETRA>                                 M_pressureExporter;
    std::string                                       M_name;
    shp<BlockVector>                                  M_extrapolatedSolution;
    shp<BlockMatrix>                                  M_mass;
    shp<BlockMatrix>                                  M_massPressure;
    shp<BlockMatrix>                                  M_stiffness;
    shp<BlockMatrix>                                  M_divergence;
    shp<FESPACE>                                      M_velocityFESpace;
    shp<FESPACE>                                      M_pressureFESpace;
    shp<FESPACE>                                      M_displacementFESpace;
    shp<ETFESPACE3>                                   M_velocityFESpaceETA;
    shp<ETFESPACE3>                                   M_displacementFESpaceETA;
    shp<ETFESPACE1>                                   M_pressureFESpaceETA;
    double                                            M_density;
    double                                            M_viscosity;
    shp<MATRIXEPETRA>                                 M_massWall;

    // first index is face flag
    std::map<unsigned int, shp<VECTOREPETRA>>         M_flowRateVectors;
    std::map<unsigned int, shp<BlockMatrix>>          M_flowRateJacobians;
    std::map<unsigned int, shp<BlockMatrix>>          M_additionalOutletMatrices;

    std::string                                       M_velocityOrder;
    std::string                                       M_pressureOrder;
    std::string                                       M_displacementOrder;
    shp<RBBases>                                      M_bases;

    shp<VECTOREPETRA>                                 M_xs;
    shp<VECTOREPETRA>                                 M_ys;
    shp<VECTOREPETRA>                                 M_zs;

    bool                                              M_addNoSlipBC;
};

}

#endif // STOKESASSEMBLERFE_HPP
