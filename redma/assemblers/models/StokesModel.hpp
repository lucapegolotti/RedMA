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

#ifndef STOKESMODEL_HPP
#define STOKESMODEL_HPP

#include <redma/RedMA.hpp>
#include <redma/RedMA.hpp>
#include <redma/assemblers/abstract/aAssembler.hpp>
#include <redma/array/DistributedVector.hpp>
#include <redma/array/SparseMatrix.hpp>
#include <redma/geometry/BuildingBlock.hpp>

// #include <redma/reduced_basis/BlockMDEIM.hpp>

#include <lifev/eta/expression/Integrate.hpp>
#include <lifev/core/filter/Exporter.hpp>
#include <lifev/core/filter/ExporterVTK.hpp>
#include <lifev/core/filter/ExporterHDF5.hpp>

namespace RedMA
{

class StokesModel
{
public:
    StokesModel(const DataContainer& data, shp<TreeNode> treeNode);

    shp<aVector> getForcingTerm(const double& time) const;

    virtual shp<aMatrix> assembleReducedStiffness(shp<BCManager> bcManager);

    virtual shp<aMatrix> assembleReducedMass(shp<BCManager> bcManager);

    virtual shp<aMatrix> assembleReducedDivergence(shp<BCManager> bcManager);

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

    void exportNorms(const double& t, shp<VECTOREPETRA> velocity, shp<VECTOREPETRA> pressure);

    void initializePythonStructures();

    void computeWallShearStress(shp<VECTOREPETRA> velocity,
                                shp<VECTOREPETRA> WSS,
                                EPETRACOMM comm);

    void applyDirichletBCsMatrix(shp<BCManager> bcManager,
                                  shp<aMatrix> matrix, double diagCoeff);

    void initializeVelocityFESpace(EPETRACOMM comm);

    void initializePressureFESpace(EPETRACOMM comm);

    void setVelocityOrder(std::string velocityOrder) {M_velocityOrder = velocityOrder;}

    void setPressureOrder(std::string pressureOrder) {M_pressureOrder = pressureOrder;}

protected:

    shp<BlockVector> buildZeroVector() const;

    std::string                                       M_name;
    DataContainer                                     M_dataContainer;
    shp<TreeNode>                                     M_treeNode;
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

    // rb structures
    // shp<BlockMDEIM>                                   M_mdeimMass;
    // shp<BlockMDEIM>                                   M_mdeimStiffness;
    // shp<BlockMDEIM>                                   M_mdeimDivergence;
    shp<RBBases>                                      M_bases;
    // PyObject*                                      M_pFunc;
    // PyObject*                                      M_pModule;

    shp<VECTOREPETRA>                                 M_xs;
    shp<VECTOREPETRA>                                 M_ys;
    shp<VECTOREPETRA>                                 M_zs;

    bool                                              M_addNoSlipBC;
};

}

#endif // STOKESMODEL_HPP
