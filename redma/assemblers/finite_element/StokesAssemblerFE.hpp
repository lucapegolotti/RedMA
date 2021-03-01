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
// #include <redma/assemblers/models/StokesModel.hpp>

namespace RedMA
{

/*! \brief Finite element assembler of the Stokes problem.
 *
 */
class StokesAssemblerFE : public aAssemblerFE
{
public:
    /*! \brief Constructor taking a datafile and a TreeNode as argument.
     *
     * \param datafile The datafile.
     * \param datafile The TreeNode encoding the physical domain.
     */
    StokesAssemblerFE(const DataContainer& data, shp<TreeNode> treeNode);

    /// Virtual setup function.
    virtual void setup() override;

    /*! Virtual export solution.
     *
     * \param time Current time.
     * \param sol Current solution.
     */
    virtual void exportSolution(const double& time,
                                const shp<aVector>& sol) override;

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

    /*! \brief Virtual getter for mass matrix jacobian.
     *
     * \param time Current time.
     * \param sol Current solution.
     * \return Shared pointer to aMatrix of the mass matrix jacobian.
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
     * \return Shared pointer to aMatrix of the right-hand side jacobian.
     */
    virtual shp<aMatrix> getJacobianRightHandSide(const double& time,
                                                  const shp<aVector>& sol) override;

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
     * \return Shared pointer to desired finite element space.
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

    virtual inline shp<ETFESPACE3> getETFESpaceCoupling() const override
    {
        return this->M_velocityFESpaceETA;
    }

    // virtual inline shp<ETFESPACE1> getETFESpaceSecondary() const override
    // {
    //     return this->M_pressureFESpaceETA;
    // }

    void applyDirichletBCsMatrix(shp<aMatrix> matrix, double diagCoeff) const override;

    void apply0DirichletBCs(shp<aVector> vector) const override;

    void applyDirichletBCs(const double& time, shp<aVector> vector) const override;

    virtual shp<FESPACE> getFEspace(unsigned int index) const override;

    shp<ETFESPACE3> getVelocityETFEspace() const {return M_velocityFESpaceETA;}

    virtual std::vector<shp<aMatrix>> getMatrices() const override;

    virtual shp<aMatrix> assembleMatrix(const unsigned int& index) override;

    virtual shp<aMatrix> getNorm(const unsigned int& fieldIndex, bool bcs = true) override;

    virtual shp<aMatrix> getConstraintMatrix() override;

    // virtual void setMDEIMs(shp<MDEIMManager> mdeimManager) override;

    void setExtrapolatedSolution(const shp<aVector>& exSol) override;

    virtual void applyPiola(shp<aVector> solution, bool inverse) override;

    void addNeumannBCs(double time, shp<aVector> sol, shp<aVector> rhs);

    void setVelocityOrder(std::string velocityOrder) {M_velocityOrder = velocityOrder;}

    void setPressureOrder(std::string pressureOrder) {M_pressureOrder = pressureOrder;}

    // void addNeumannBCs(shp<aVector>& input, const double& time,
    //                    const shp<aVector>& sol);

    shp<aMatrix> assembleStiffness(shp<BCManager> bcManager);

    shp<aMatrix> assembleMass(shp<BCManager> bcManager);

    shp<aMatrix> assembleDivergence(shp<BCManager> bcManager);

    std::map<unsigned int, double> computeFlowRates(shp<aVector> sol,
                                                    bool verbose = false);

    // these are int_{Gamma_N} phi_i * n inlets (just to check that it is preserved)
    // and outlets (for boundary conditions)
    void assembleFlowRateVectors();

    // compute the jacobian of int_{Gamma_N}(int_{Gamma_N} u_i * n)phi_j
    void assembleFlowRateJacobians(shp<BCManager> bcManager);

    shp<VECTOREPETRA> assembleFlowRateVector(const GeometricFace& face);

    shp<MATRIXEPETRA> assembleFlowRateJacobian(const GeometricFace& face);

    void addBackFlowStabilization(shp<aVector>& input,
                                  shp<aVector> sol,
                                  const unsigned int& faceFlag);

    void exportNorms(double t, shp<VECTOREPETRA> velocity, shp<VECTOREPETRA> pressure);

    void initializePythonStructures();

    void computeWallShearStress(shp<VECTOREPETRA> velocity,
                                shp<VECTOREPETRA> WSS,
                                EPETRACOMM comm);

    void applyDirichletBCsMatrix(shp<BCManager> bcManager,
                                  shp<aMatrix> matrix, double diagCoeff);

    void initializeVelocityFESpace(EPETRACOMM comm);

    void initializePressureFESpace(EPETRACOMM comm);

    shp<aVector> getForcingTerm(const double& time) const;

    inline double getDensity() {return M_density;}

    inline double getViscosity() {return M_viscosity;}

protected:
    shp<BlockVector> buildZeroVector() const;

    shp<LifeV::Exporter<MESH>>                        M_exporter;
    shp<VECTOREPETRA>                                 M_velocityExporter;
    shp<VECTOREPETRA>                                 M_WSSExporter;
    shp<VECTOREPETRA>                                 M_pressureExporter;
    std::string                                       M_name;
    shp<BlockVector>                                  M_extrapolatedSolution;
    shp<BlockMatrix>                                  M_mass;
    shp<BlockMatrix>                                  M_stiffness;
    shp<BlockMatrix>                                  M_divergence;
    shp<FESPACE>                                      M_velocityFESpace;
    shp<FESPACE>                                      M_pressureFESpace;
    shp<ETFESPACE3>                                   M_velocityFESpaceETA;
    shp<ETFESPACE1>                                   M_pressureFESpaceETA;
    double                                            M_density;
    double                                            M_viscosity;
    shp<MATRIXEPETRA>                                 M_massWall;
    shp<SparseMatrix>                                 M_massVelocity;
    shp<SparseMatrix>                                 M_massPressure;
    // first index is face flag
    std::map<unsigned int, shp<VECTOREPETRA>>         M_flowRateVectors;
    std::map<unsigned int, shp<BlockMatrix>>          M_flowRateJacobians;

    std::string                                       M_velocityOrder;
    std::string                                       M_pressureOrder;
    shp<RBBases>                                      M_bases;

    shp<VECTOREPETRA>                                 M_xs;
    shp<VECTOREPETRA>                                 M_ys;
    shp<VECTOREPETRA>                                 M_zs;
};

}

#endif // STOKESASSEMBLERFE_HPP
