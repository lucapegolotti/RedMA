// // Reduced Modeling of Arteries (RedMA)
// // Copyright (C) 2019  Luca Pegolotti
// //
// // RedMA is free software: you can redistribute it and/or modify
// // it under the terms of the GNU General Public License as published by
// // the Free Software Foundation, either version 3 of the License, or
// // (at your option) any later version.
// //
// // RedMA is distributed in the hope that it will be useful,
// // but WITHOUT ANY WARRANTY; without even the implied warranty of
// // MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// // GNU General Public License for more details.
// //
// // You should have received a copy of the GNU General Public License
// // along with this program.  If not, see <http://www.gnu.org/licenses/>.
//
// #ifndef STOKESASSEMBLER_HPP
// #define STOKESASSEMBLER_HPP
//
// #include <redma/assemblers/abstract/aAssembler.hpp>
// #include <redma/array/DistributedVector.hpp>
// #include <redma/array/SparseMatrix.hpp>
// #include <redma/geometry/BuildingBlock.hpp>
//
// #include <redma/reduced_basis/BlockMDEIM.hpp>
//
// #include <lifev/eta/expression/Integrate.hpp>
// #include <lifev/core/filter/Exporter.hpp>
// #include <lifev/core/filter/ExporterVTK.hpp>
// #include <lifev/core/filter/ExporterHDF5.hpp>
//
// #include <boost/filesystem.hpp>
// // #include <Python.h>
//
// // class BlockMDEIM;
//
// namespace RedMA
// {
//
// class StokesAssembler : public aAssembler
// {
// public:
//     StokesAssembler(const DataContainer& data, SHP(TreeNode) treeNode);
//
//     virtual void setup() override;
//
//     virtual void exportSolution(const double& t,
//                                 const BlockVector& sol) override;
//
//     virtual void postProcess(const double& t, const BlockVector& sol) override;
//
//     virtual BlockMatrix getMass(const double& time,
//                                 const BlockVector& sol) override;
//
//     virtual BlockMatrix getMassJacobian(const double& time,
//                                         const BlockVector& sol) override;
//
//     virtual BlockVector getRightHandSide(const double& time,
//                                       const BlockVector& sol) override;
//
//     virtual BlockMatrix getJacobianRightHandSide(const double& time,
//                                       const BlockVector& sol) override;
//
//     virtual BlockVector getZeroVector() const override;
//
//     virtual BlockVector getLifting(const double& time) const override;
//
//     BlockVector getForcingTerm(const double& time) const;
//
//     void addNeumannBCs(BlockVector& input, const double& time,
//                        const BlockVector& sol);
//
//     void initializeFEspaces() override;
//
//     BlockMatrix assembleStiffness(BlockMDEIMStructure* structure = nullptr);
//
//     BlockMatrix assembleReducedStiffness(BlockMDEIMStructure* structure = nullptr);
//
//     BlockMatrix assembleMass(BlockMDEIMStructure* structure = nullptr);
//
//     BlockMatrix assembleReducedMass(BlockMDEIMStructure* structure = nullptr);
//
//     BlockMatrix assembleDivergence(BlockMDEIMStructure* structure = nullptr);
//
//     BlockMatrix assembleReducedDivergence(BlockMDEIMStructure* structure = nullptr);
//
//     BlockVector getFELifting(const double& time) const override;
//
//     std::map<unsigned int, double> computeFlowRates(const BlockVector& sol,
//                                                     bool verbose = false);
//
//     // these are int_{Gamma_N} phi_i * n inlets (just to check that it is preserved)
//     // and outlets (for boundary conditions)
//     void assembleFlowRateVectors();
//
//     // compute the jacobian of int_{Gamma_N}(int_{Gamma_N} u_i * n)phi_j
//     void assembleFlowRateJacobians();
//
//     SHP(VECTOREPETRA) assembleFlowRateVector(const GeometricFace& face);
//
//     SHP(MATRIXEPETRA) assembleFlowRateJacobian(const GeometricFace& face);
//
//     void setExporter() override;
//
//     virtual inline SHP(FESPACE) getFESpaceBCs() const override
//     {
//         return M_velocityFESpace;
//     }
//
//     virtual inline unsigned int getComponentBCs() const override {return 0;}
//
//     virtual inline SHP(ETFESPACE3) getETFESpaceCoupling() const override
//     {
//         return M_velocityFESpaceETA;
//     }
//
//     virtual inline SHP(ETFESPACE1) getETFESpaceSecondary() const override
//     {
//         return M_pressureFESpaceETA;
//     }
//
//     void addBackFlowStabilization(BlockVector input,
//                                   const BlockVector& sol,
//                                   const unsigned int& faceFlag);
//
//     void apply0DirichletBCsMatrix(BlockMatrix& matrix, double diagCoeff) const override;
//
//     void apply0DirichletBCs(BlockVector& vector) const override;
//
//     void applyDirichletBCs(const double& time, BlockVector& vector) const override;
//
//     virtual inline SHP(FESPACE) getFEspace(unsigned int index) const override;
//
//     virtual std::vector<BlockMatrix> getMatrices() const override;
//
//     virtual BlockMatrix assembleMatrix(const unsigned int& index,
//                                        BlockMDEIMStructure* structure = nullptr) override;
//
//     virtual SparseMatrix getNorm(const unsigned int& fieldIndex, bool bcs = true) override;
//
//     virtual std::shared_ptr<aMatrix> getConstraintMatrix() override;
//
//     virtual void setMDEIMs(SHP(MDEIMManager) mdeimManager) override;
//
//     virtual void setRBBases(SHP(RBBasesManager) rbManager) override;
//
//     virtual SHP(RBBases) getRBBases() const override {return M_bases;}
//
//     virtual BlockVector convertFunctionRBtoFEM(BlockVector rbSolution) const override;
//
//     void exportNorms(double t);
//
//     virtual void RBsetup() override;
//
//     void setExtrapolatedSolution(const BlockVector& exSol) override;
//
//     void initializePythonStructures();
//
//     virtual void applyPiola(BlockVector solution, bool inverse) override;
//
//     void computeWallShearStress(SHP(VECTOREPETRA) velocity);
//
// protected:
//
//     BlockMatrix                                       M_mass;
//     BlockMatrix                                       M_stiffness;
//     BlockMatrix                                       M_divergence;
//     SHP(FESPACE)                                      M_velocityFESpace;
//     SHP(FESPACE)                                      M_pressureFESpace;
//     SHP(ETFESPACE3)                                   M_velocityFESpaceETA;
//     SHP(ETFESPACE1)                                   M_pressureFESpaceETA;
//     double                                            M_density;
//     double                                            M_viscosity;
//     SHP(VECTOREPETRA)                                 M_velocityExporter;
//     SHP(VECTOREPETRA)                                 M_WSSExporter;
//     SHP(VECTOREPETRA)                                 M_pressureExporter;
//     SHP(MATRIXEPETRA)                                 M_massWall;
//     MatrixEp                                          M_massVelocity;
//     MatrixEp                                          M_massPressure;
//     SHP(LifeV::Exporter<MESH>)                        M_exporter;
//     // first index is face flag
//     std::map<unsigned int, SHP(VECTOREPETRA)>         M_flowRateVectors;
//     std::map<unsigned int, BlockMatrix>               M_flowRateJacobians;
//
//     BlockVector                                       M_extrapolatedSolution;
//
//     // rb structures
//     SHP(BlockMDEIM)                                   M_mdeimMass;
//     SHP(BlockMDEIM)                                   M_mdeimStiffness;
//     SHP(BlockMDEIM)                                   M_mdeimDivergence;
//     SHP(RBBases)                                      M_bases;
//     // PyObject*                                      M_pFunc;
//     // PyObject*                                      M_pModule;
//
//     SHP(VECTOREPETRA)                                 M_xs;
//     SHP(VECTOREPETRA)                                 M_ys;
//     SHP(VECTOREPETRA)                                 M_zs;
// };
//
// }
//
// #include "StokesAssembler_imp.hpp"
//
// #endif // STOKESASSEMBLER_HPP
