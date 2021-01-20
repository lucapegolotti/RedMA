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

#ifndef STOKESASSEMBLERRB_HPP
#define STOKESASSEMBLERRB_HPP

#include <redma/assemblers/abstract/aAssemblerFE.hpp>
#include <redma/assemblers/finite_element/StokesAssemblerFE.hpp>
#include <redma/assemblers/abstract/aAssemblerRB.hpp>
#include <redma/reduced_basis/RBBases.hpp>

namespace RedMA
{

class StokesAssemblerRB : public aAssemblerRB, public StokesModel
{
public:
    StokesAssemblerRB(const DataContainer& data, shp<TreeNode> treeNode);

    virtual void setup() override;

    virtual void exportSolution(const double& t,
                                const shp<aVector>& sol) override;

    virtual void postProcess(const double& t, const shp<aVector>& sol) override;

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

    // virtual inline shp<ETFESPACE1> getETFESpaceSecondary() const override
    // {
    //     return this->M_pressureFESpaceETA;
    // }

    void applyDirichletBCsMatrix(shp<aMatrix> matrix, double diagCoeff) const override;

    void apply0DirichletBCs(shp<aVector> vector) const override;

    void applyDirichletBCs(const double& time, shp<aVector> vector) const override;

    virtual inline shp<FESPACE> getFEspace(unsigned int index) const override {}

    virtual std::vector<shp<aMatrix>> getMatrices() const override;

    virtual shp<aMatrix> assembleMatrix(const unsigned int& index) override;

    // virtual void setMDEIMs(shp<MDEIMManager> mdeimManager) override {throw new Exception("setMDEIMs method not implemented for RB");}

    void setExtrapolatedSolution(const shp<aVector>& exSol) override {throw new Exception("setExtrapolatedSolution method not implemented for RB");}

    virtual void applyPiola(shp<aVector> solution, bool inverse) override;

    virtual void RBsetup() override;

    virtual shp<RBBases> getRBBases() const override;

    virtual void setRBBases(shp<RBBasesManager> rbManager) override;

    virtual shp<aVector> convertFunctionRBtoFEM(shp<aVector> rbSolution) const override;

protected:
    shp<LifeV::Exporter<MESH>>                        M_exporter;
    shp<VECTOREPETRA>                                 M_velocityExporter;
    shp<VECTOREPETRA>                                 M_WSSExporter;
    shp<VECTOREPETRA>                                 M_pressureExporter;
    std::string                                       M_name;
    shp<BlockVector>                                  M_extrapolatedSolution;
    shp<RBBases>                                      M_bases;
    shp<BlockMatrix>                                  M_reducedMass;
    shp<BlockMatrix>                                  M_reducedDivergence;
    shp<BlockMatrix>                                  M_reducedStiffness;
};

}

#endif // STOKESASSEMBLERRB_HPP
