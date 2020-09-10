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

#ifndef STOKESASSEMBLERFEM_HPP
#define STOKESASSEMBLERFEM_HPP

#include <redma/assemblers/abstract/aAssemblerRB.hpp>
#include <redma/assemblers/models/StokesAssembler.hpp>

namespace RedMA
{

class StokesAssemblerFEM : public aAssemblerFE, public StokesAssembler
{
public:
    StokesAssembler(const DataContainer& data, SHP(TreeNode) treeNode);

    virtual void setup() override;

    virtual void exportSolution(const double& t,
                                const BlockVector& sol) override;

    virtual void postProcess(const double& t, const BlockVector& sol) override;

    virtual BlockMatrix getMass(const double& time,
                                const BlockVector& sol) override;

    virtual BlockMatrix getMassJacobian(const double& time,
                                        const BlockVector& sol) override;

    virtual BlockVector getRightHandSide(const double& time,
                                      const BlockVector& sol) override;

    virtual BlockMatrix getJacobianRightHandSide(const double& time,
                                      const BlockVector& sol) override;

    virtual BlockVector getZeroVector() const override;

    virtual BlockVector getLifting(const double& time) const override;

    void initializeFEspaces() override;

    BlockVector getFELifting(const double& time) const override

    void setExporter() override;

    virtual inline SHP(FESPACE) getFESpaceBCs() const override
    {
        return this->M_velocityFESpace;
    }

    virtual inline unsigned int getComponentBCs() const override {return 0;}

    virtual inline SHP(ETFESPACE3) getETFESpaceCoupling() const override
    {
        return this->M_velocityFESpaceETA;
    }

    virtual inline SHP(ETFESPACE1) getETFESpaceSecondary() const override
    {
        return this->M_pressureFESpaceETA;
    }

    void apply0DirichletBCsMatrix(BlockMatrix& matrix, double diagCoeff) const override;

    void apply0DirichletBCs(BlockVector& vector) const override;

    void applyDirichletBCs(const double& time, BlockVector& vector) const override;

    virtual inline SHP(FESPACE) getFEspace(unsigned int index) const override;

    virtual std::vector<BlockMatrix> getMatrices() const override;

    virtual BlockMatrix assembleMatrix(const unsigned int& index,
                                       BlockMDEIMStructure* structure = nullptr) override;

    virtual SparseMatrix getNorm(const unsigned int& fieldIndex, bool bcs = true) override;

    virtual std::shared_ptr<aMatrix> getConstraintMatrix() override;

    virtual void setMDEIMs(SHP(MDEIMManager) mdeimManager) override;

    void setExtrapolatedSolution(const BlockVector& exSol) override;

    virtual void applyPiola(BlockVector solution, bool inverse) override;

protected:

    SHP(LifeV::Exporter<MESH>)                        M_exporter;

};

}

#endif // STOKESASSEMBLERFEM_HPP
