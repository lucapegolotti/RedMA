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
#include <redma/assemblers/models/StokesModel.hpp>

namespace RedMA
{

class StokesAssemblerFE : public aAssemblerFE, public StokesModel
{
public:
    StokesAssemblerFE(const DataContainer& data, shp<TreeNode> treeNode);

    virtual void setup() override;

    virtual void exportSolution(const double& t,
                                const shp<aVector>& sol) override;

    virtual void postProcess(const double& t, const double &dt, const shp<aVector>& sol) override;

    virtual shp<aMatrix> getMass(const double& time,
                                 const shp<aVector>& sol) override;

    virtual shp<aMatrix> getMassJacobian(const double& time,
                                         const shp<aVector>& sol) override;

    virtual shp<aVector> getRightHandSide(const double& time,
                                              const shp<aVector>& sol) override;

    virtual shp<aMatrix> getJacobianRightHandSide(const double& time,
                                                      const shp<aVector>& sol) override;

    virtual shp<aVector> getZeroVector() const override;

    virtual shp<aVector> getLifting(const double& time) const override;

    virtual shp<aVector> getFELifting(const double& time) const override;

    void initializeFEspaces() override;

    void setExporter() override;

    virtual inline shp<FESPACE> getFESpaceBCs() const override
    {
        return this->M_velocityFESpace;
    }

    virtual inline unsigned int getComponentBCs() const override {return 0;}

    virtual inline shp<ETFESPACE3> getETFESpaceCoupling() const override
    {
        return this->M_velocityFESpaceETA;
    }

    void applyDirichletBCsMatrix(shp<aMatrix> matrix, double diagCoeff) const override;

    void apply0DirichletBCs(shp<aVector> vector) const override;

    void applyDirichletBCs(const double& time, shp<aVector> vector) const override;

    virtual shp<FESPACE> getFEspace(unsigned int index) const override;

    virtual std::vector<shp<aMatrix>> getMatrices() const override;

    virtual shp<aMatrix> assembleMatrix(const unsigned int& index) override;

    virtual shp<aMatrix> getNorm(const unsigned int& fieldIndex, bool bcs = true) override;

    virtual shp<aMatrix> getConstraintMatrix() override;

    // virtual void setMDEIMs(shp<MDEIMManager> mdeimManager) override;

    void setExtrapolatedSolution(const shp<aVector>& exSol) override;

    virtual void applyPiola(shp<aVector> solution, bool inverse) override;

    void addNeumannBCs(double time, shp<aVector> sol, shp<aVector> rhs);

protected:
    shp<LifeV::Exporter<MESH>>                        M_exporter;
    shp<VECTOREPETRA>                                 M_velocityExporter;
    shp<VECTOREPETRA>                                 M_WSSExporter;
    shp<VECTOREPETRA>                                 M_pressureExporter;
    std::string                                       M_name;
    shp<BlockVector>                                  M_extrapolatedSolution;
};

}

#endif // STOKESASSEMBLERFE_HPP
