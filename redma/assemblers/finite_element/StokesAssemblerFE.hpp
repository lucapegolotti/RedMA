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

#include <redma/assemblers/abstract/aAssemblerFE.hpp>
#include <redma/assemblers/models/StokesAssembler.hpp>

namespace RedMA
{

class StokesAssemblerFE : public aAssemblerFE, public StokesAssembler
{
public:
    StokesAssemblerFE(const DataContainer& data, SHP(TreeNode) treeNode);

    virtual void setup() override;

    virtual void exportSolution(const double& t,
                                const SHP(aVector)& sol) override;

    virtual void postProcess(const double& t, const SHP(aVector)& sol) override;

    virtual SHP(aMatrix) getMass(const double& time,
                                     const SHP(aVector)& sol) override;

    virtual SHP(aMatrix) getMassJacobian(const double& time,
                                             const SHP(aVector)& sol) override;

    virtual SHP(aVector) getRightHandSide(const double& time,
                                              const SHP(aVector)& sol) override;

    virtual SHP(aMatrix) getJacobianRightHandSide(const double& time,
                                                      const SHP(aVector)& sol) override;

    virtual SHP(aVector) getZeroVector() const override;

    virtual SHP(aVector) getLifting(const double& time) const override;

    virtual SHP(aVector) getFELifting(const double& time) const override;

    void initializeFEspaces() override;

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

    void apply0DirichletBCsMatrix(SHP(aMatrix) matrix, double diagCoeff) const override;

    void apply0DirichletBCs(SHP(aVector) vector) const override;

    void applyDirichletBCs(const double& time, SHP(aVector) vector) const override;

    virtual inline SHP(FESPACE) getFEspace(unsigned int index) const override;

    virtual std::vector<SHP(aMatrix)> getMatrices() const override;

    virtual SHP(aMatrix) assembleMatrix(const unsigned int& index,
                                            BlockMDEIMStructure* structure = nullptr) override;

    virtual SHP(aMatrix) getNorm(const unsigned int& fieldIndex, bool bcs = true) override;

    virtual SHP(aMatrix) getConstraintMatrix() override;

    virtual void setMDEIMs(SHP(MDEIMManager) mdeimManager) override;

    void setExtrapolatedSolution(const SHP(aVector)& exSol) override;

    virtual void applyPiola(SHP(aVector) solution, bool inverse) override;

protected:
    SHP(LifeV::Exporter<MESH>)                        M_exporter;
    SHP(VECTOREPETRA)                                 M_velocityExporter;
    SHP(VECTOREPETRA)                                 M_WSSExporter;
    SHP(VECTOREPETRA)                                 M_pressureExporter;
    std::string                                       M_name;
    SHP(BlockVector)                                  M_extrapolatedSolution;
};

}

#endif // STOKESASSEMBLERFE_HPP
