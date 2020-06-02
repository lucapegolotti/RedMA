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

#ifndef STOKESASSEMBLER_HPP
#define STOKESASSEMBLER_HPP

#include <redma/solver/assemblers/aAssembler.hpp>
#include <redma/solver/array/VectorEp.hpp>
#include <redma/solver/array/MatrixEp.hpp>
#include <redma/geometry/BuildingBlock.hpp>

#include <redma/reduced_basis/BlockMDEIM.hpp>

#include <lifev/eta/expression/Integrate.hpp>
#include <lifev/core/filter/Exporter.hpp>
#include <lifev/core/filter/ExporterVTK.hpp>
#include <lifev/core/filter/ExporterHDF5.hpp>

#include <boost/filesystem.hpp>
#include <Python.h>

// class BlockMDEIM;

namespace RedMA
{

template <class InVectorType, class InMatrixType>
class StokesAssembler : public aAssembler<InVectorType, InMatrixType>
{
public:
    StokesAssembler(const DataContainer& data, SHP(TreeNode) treeNode);

    virtual void setup() override;

    virtual void exportSolution(const double& t,
                                const BlockVector<InVectorType>& sol) override;

    virtual void postProcess(const double& t, const BlockVector<InVectorType>& sol) override;

    virtual BlockMatrix<InMatrixType> getMass(const double& time,
                                      const BlockVector<InVectorType>& sol) override;

    virtual BlockMatrix<InMatrixType> getMassJacobian(const double& time,
                                      const BlockVector<InVectorType>& sol) override;

    virtual BlockVector<InVectorType> getRightHandSide(const double& time,
                                      const BlockVector<InVectorType>& sol) override;

    virtual BlockMatrix<InMatrixType> getJacobianRightHandSide(const double& time,
                                      const BlockVector<InVectorType>& sol) override;

    virtual BlockVector<InVectorType> getZeroVector() const override;

    virtual BlockVector<InVectorType> getLifting(const double& time) const override;

    BlockVector<InVectorType> getForcingTerm(const double& time) const;

    void addNeumannBCs(BlockVector<FEVECTOR>& input, const double& time,
                       const BlockVector<InVectorType>& sol);

    void initializeFEspaces();

    BlockMatrix<InMatrixType> assembleStiffness(BlockMDEIMStructure* structure = nullptr);

    BlockMatrix<MatrixEp> assembleReducedStiffness(BlockMDEIMStructure* structure = nullptr);

    BlockMatrix<InMatrixType> assembleMass(BlockMDEIMStructure* structure = nullptr);

    BlockMatrix<MatrixEp> assembleReducedMass(BlockMDEIMStructure* structure = nullptr);

    BlockMatrix<InMatrixType> assembleDivergence(BlockMDEIMStructure* structure = nullptr);

    BlockMatrix<MatrixEp> assembleReducedDivergence(BlockMDEIMStructure* structure = nullptr);

    BlockVector<FEVECTOR> getFELifting(const double& time) const override;

    std::map<unsigned int, double> computeFlowRates(const BlockVector<InVectorType>& sol,
                                                    bool verbose = false);

    // these are int_{Gamma_N} phi_i * n inlets (just to check that it is preserved)
    // and outlets (for boundary conditions)
    void assembleFlowRateVectors();

    // compute the jacobian of int_{Gamma_N}(int_{Gamma_N} u_i * n)phi_j
    void assembleFlowRateJacobians();

    SHP(VECTOREPETRA) assembleFlowRateVector(const GeometricFace& face);

    SHP(MATRIXEPETRA) assembleFlowRateJacobian(const GeometricFace& face);

    void setExporter() override;

    virtual inline SHP(FESPACE) getFESpaceBCs() const override
    {
        return M_velocityFESpace;
    }

    virtual inline unsigned int getComponentBCs() const override {return 0;}

    virtual inline SHP(ETFESPACE3) getETFESpaceCoupling() const override
    {
        return M_velocityFESpaceETA;
    }

    virtual inline SHP(ETFESPACE1) getETFESpaceSecondary() const override
    {
        return M_pressureFESpaceETA;
    }

    void addBackFlowStabilization(BlockVector<InVectorType> input,
                                  const BlockVector<InVectorType>& sol,
                                  const unsigned int& faceFlag);

    void apply0DirichletBCsMatrix(BlockMatrix<InMatrixType>& matrix, double diagCoeff) const override;

    void apply0DirichletBCs(BlockVector<InVectorType>& vector) const override;

    void applyDirichletBCs(const double& time, BlockVector<InVectorType>& vector) const override;

    virtual inline SHP(FESPACE) getFEspace(unsigned int index) const override;

    virtual std::vector<BlockMatrix<InMatrixType>> getMatrices() const override;

    virtual BlockMatrix<InMatrixType> assembleMatrix(const unsigned int& index,
                                                     BlockMDEIMStructure* structure = nullptr) override;

    virtual MatrixEp getNorm(const unsigned int& fieldIndex, bool bcs = true) override;

    virtual InMatrixType getConstraintMatrix() override;

    virtual void setMDEIMs(SHP(MDEIMManager) mdeimManager) override;

    virtual void setRBBases(SHP(RBBasesManager) rbManager) override;

    virtual SHP(RBBases) getRBBases() const override {return M_bases;}

    virtual BlockVector<FEVECTOR> convertFunctionRBtoFEM(BlockVector<RBVECTOR> rbSolution) const override;

    void exportNorms(double t);

    virtual void RBsetup() override;

    void setExtrapolatedSolution(const BlockVector<InVectorType>& exSol) override;

    void initializePythonStructures();

protected:

    BlockMatrix<InMatrixType>                                       M_mass;
    BlockMatrix<InMatrixType>                                       M_stiffness;
    BlockMatrix<InMatrixType>                                       M_divergence;
    SHP(FESPACE)                                                    M_velocityFESpace;
    SHP(FESPACE)                                                    M_pressureFESpace;
    SHP(ETFESPACE3)                                                 M_velocityFESpaceETA;
    SHP(ETFESPACE1)                                                 M_pressureFESpaceETA;
    double                                                          M_density;
    double                                                          M_viscosity;
    SHP(VECTOREPETRA)                                               M_velocityExporter;
    SHP(VECTOREPETRA)                                               M_pressureExporter;
    SHP(VECTOREPETRA)                                               M_WSSExporter;
    MatrixEp                                                        M_massVelocity;
    MatrixEp                                                        M_massPressure;
    SHP(LifeV::Exporter<MESH>)                                      M_exporter;
    // first index is face flag
    std::map<unsigned int, SHP(VECTOREPETRA)>                       M_flowRateVectors;
    std::map<unsigned int, BlockMatrix<InMatrixType>>               M_flowRateJacobians;

    BlockVector<InVectorType>                                       M_extrapolatedSolution;

    // rb structures
    SHP(BlockMDEIM)                                                 M_mdeimMass;
    SHP(BlockMDEIM)                                                 M_mdeimStiffness;
    SHP(BlockMDEIM)                                                 M_mdeimDivergence;
    SHP(RBBases)                                                    M_bases;
    PyObject*                                                       M_pFunc;
    PyObject*                                                       M_pModule;
};

}

#include "StokesAssembler_imp.hpp"

#endif // STOKESASSEMBLER_HPP
