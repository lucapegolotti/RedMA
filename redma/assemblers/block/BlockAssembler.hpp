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

#ifndef BLOCKASSEMBLER_HPP
#define BLOCKASSEMBLER_HPP

#include <redma/RedMA.hpp>
#include <redma/assemblers/abstract/aAssembler.hpp>
#include <redma/assemblers/abstract/aAssemblerFE.hpp>
#include <redma/assemblers/abstract/aAssemblerRB.hpp>
#include <redma/assemblers/AssemblerFactory.hpp>
#include <redma/assemblers/coupling/InterfaceAssembler.hpp>
#include <redma/assemblers/coupling/InletInflowAssembler.hpp>

#include <redma/reduced_basis/RBBasesManager.hpp>

#include <redma/geometry/TreeStructure.hpp>

namespace RedMA
{

class BlockAssembler : public aAssembler
{
    typedef DefaultAssemblersLibrary  DefaultAssemblers;

public:
    BlockAssembler() {}

    BlockAssembler(const DataContainer& data, const TreeStructure& tree,
                   shp<DefaultAssemblers> defAssemblers = nullptr);

    // virtual ~BlockAssembler() {}

    virtual void setup() override;

    virtual void exportSolution(const double& t,
                                const shp<aVector>& sol) override;

    virtual void postProcess(const double& t,
                             const shp<aVector>& sol) override;

    virtual shp<aMatrix> getMass(const double& time,
                                const shp<aVector>& sol) override;

    virtual shp<aMatrix> getMassJacobian(const double& time,
                                        const shp<aVector>& sol) override;

    virtual shp<aVector> getRightHandSide(const double& time,
                                         const shp<aVector>& sol) override;

    virtual shp<aMatrix> getJacobianRightHandSide(const double& time,
                                                 const shp<aVector>& sol) override;

    virtual shp<aVector> getLifting(const double& time) const override;

    virtual shp<aVector> getZeroVector() const override;

    virtual void applyDirichletBCsMatrix(shp<aMatrix> matrix, double diagCoeff) const override;

    virtual void apply0DirichletBCs(shp<aVector> vector) const override;

    virtual void setExporter() override;

    virtual void applyDirichletBCs(const double& time, shp<aVector> vector) const override;

    virtual void checkStabTerm(const shp<aVector>& sol) const;

    std::map<unsigned int, std::string> getIDMeshTypeMap() const;

    inline shp<aAssembler> block(const unsigned int& index) {return M_primalAssemblers[index];}

    std::map<unsigned int, shp<aAssembler>> getAssemblersMap() const {return M_primalAssemblers;}

    std::vector<shp<InterfaceAssembler>> getDualAssemblers() const {return M_dualAssemblers;}

    shp<aVector> convertFunctionRBtoFEM(shp<aVector> rbFunction, EPETRACOMM comm) const;

    virtual void setExtrapolatedSolution(const shp<aVector>& exSol) override;

    virtual shp<aVector> getNonLinearTerm() override;

    std::map<unsigned int,std::vector<double>> getRandomizibleParametersVectors();

    virtual void applyPiola(shp<aVector> solution, bool inverse) override {}

    void applyGlobalPiola(shp<aVector> solution, bool inverse);

    virtual void initializeFEspaces() override;

    void setDefaultAssemblers(shp<DefaultAssemblers> defAssemblers) override;

protected:
    GetPot                                                        M_datafile;
    TreeStructure                                                 M_tree;
    std::map<unsigned int, shp<aAssembler>>                       M_primalAssemblers;
    std::vector<shp<InterfaceAssembler>>                          M_dualAssemblers;
    unsigned int                                                  M_numberBlocks;
    shp<MDEIMManager>                                             M_mdeimManager;
    shp<RBBasesManager>                                           M_basesManager;
};

}

#endif // BLOCKASSEMBLER_HPP
