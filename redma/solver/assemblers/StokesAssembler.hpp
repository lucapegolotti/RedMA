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

#include <lifev/eta/expression/Integrate.hpp>
#include <lifev/core/filter/Exporter.hpp>
#include <lifev/core/filter/ExporterVTK.hpp>
#include <lifev/core/filter/ExporterHDF5.hpp>

#include <boost/filesystem.hpp>

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

    virtual void postProcess(const BlockVector<InVectorType>& sol) override;

    virtual BlockMatrix<InMatrixType> getMass(const double& time,
                                      const BlockVector<InVectorType>& sol) override;

    virtual BlockVector<InVectorType> getRightHandSide(const double& time,
                                      const BlockVector<InVectorType>& sol) override;

    virtual BlockMatrix<InMatrixType> getJacobianRightHandSide(const double& time,
                                      const BlockVector<InVectorType>& sol) override;

    virtual BlockVector<InVectorType> getZeroVector() const override;

    BlockVector<FEVECTOR> computeLifting(const double& time) const;

    BlockVector<FEVECTOR> computeLiftingDt(const double& time) const;

    void addNeumannBCs(BlockVector<FEVECTOR>& input, const double& time) const;

    void initializeFEspaces();

    void assembleStiffness();

    void assembleMass();

    void assembleDivergence();

    void setExporter();

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

protected:
    void apply0DirichletBCs(BlockVector<InVectorType>& vector);

    BlockMatrix<InMatrixType>            M_mass;
    BlockMatrix<InMatrixType>            M_stiffness;
    BlockMatrix<InMatrixType>            M_divergence;
    SHP(FESPACE)                         M_velocityFESpace;
    SHP(FESPACE)                         M_pressureFESpace;
    SHP(ETFESPACE3)                      M_velocityFESpaceETA;
    SHP(ETFESPACE1)                      M_pressureFESpaceETA;
    double                               M_density;
    double                               M_viscosity;
    SHP(LifeV::VectorEpetra)             M_velocityExporter;
    SHP(LifeV::VectorEpetra)             M_pressureExporter;
    SHP(LifeV::Exporter<MESH>)           M_exporter;
};

}

#include "StokesAssembler_imp.hpp"

#endif // STOKESASSEMBLER_HPP
