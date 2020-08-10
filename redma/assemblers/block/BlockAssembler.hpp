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
    BlockAssembler(const DataContainer& data, const TreeStructure& tree,
                   SHP(DefaultAssemblers) defAssemblers = nullptr);

    virtual void setup() override;

    virtual void exportSolution(const double& t,
                                const BlockVector& sol) override;

    virtual void postProcess(const double& t,
                             const BlockVector& sol) override;

    virtual BlockMatrix getMass(const double& time,
                                const BlockVector& sol) override;

    virtual BlockMatrix getMassJacobian(const double& time,
                                        const BlockVector& sol) override;

    virtual BlockVector getRightHandSide(const double& time,
                                         const BlockVector& sol) override;

    virtual BlockMatrix getJacobianRightHandSide(const double& time,
                                                 const BlockVector& sol) override;

    virtual BlockVector getLifting(const double& time) const override;

    virtual BlockVector getZeroVector() const override;

    virtual void apply0DirichletBCsMatrix(BlockMatrix& matrix, double diagCoeff) const override;

    virtual void apply0DirichletBCs(BlockVector& vector) const override;

    virtual void setExporter() override;

    virtual void applyDirichletBCs(const double& time, BlockVector& vector) const override;

    virtual void checkStabTerm(const BlockVector& sol) const;

    std::map<unsigned int, std::string> getIDMeshTypeMap() const;

    inline SHP(aAssembler) block(const unsigned int& index) {return M_primalAssemblers[index];}

    std::map<unsigned int, SHP(aAssembler)> getAssemblersMap() const {return M_primalAssemblers;}

    std::vector<SHP(InterfaceAssembler)> getDualAssemblers() const {return M_dualAssemblers;}

    BlockVector convertFunctionRBtoFEM(BlockVector rbFunction, EPETRACOMM comm) const;

    virtual void setExtrapolatedSolution(const BlockVector& exSol) override;

    virtual BlockVector getNonLinearTerm() override;

    std::map<unsigned int,std::vector<double>> getRandomizibleParametersVectors();

    virtual void applyPiola(BlockVector solution, bool inverse) override {};

    void applyGlobalPiola(BlockVector solution, bool inverse);

    virtual void initializeFEspaces() override;

    void setDefaultAssemblers(SHP(DefaultAssemblers) defAssemblers) override;

protected:
    GetPot                                                        M_datafile;
    TreeStructure                                                 M_tree;
    std::map<unsigned int, SHP(aAssembler)>                       M_primalAssemblers;
    std::vector<SHP(InterfaceAssembler)>                          M_dualAssemblers;
    unsigned int                                                  M_numberBlocks;
    SHP(MDEIMManager)                                             M_mdeimManager;
    SHP(RBBasesManager)                                           M_basesManager;
};

}

#endif // BLOCKASSEMBLER_HPP
