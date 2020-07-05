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
#include <redma/solver/assemblers/aAssembler.hpp>
#include <redma/solver/assemblers/AssemblerFactory.hpp>
#include <redma/solver/assemblers/InterfaceAssembler.hpp>
#include <redma/solver/assemblers/InletInflowAssembler.hpp>

#include <redma/reduced_basis/RBBasesManager.hpp>

#include <redma/geometry/TreeStructure.hpp>

namespace RedMA
{

template <class InVectorType, class InMatrixType>
class BlockAssembler : public aAssembler<InVectorType, InMatrixType>
{
    typedef typename InVectorType::InnerType              VInner;
    typedef typename InMatrixType::InnerType              MInner;

public:
    BlockAssembler(const DataContainer& data, const TreeStructure& tree);

    virtual void setup() override;

    virtual void exportSolution(const double& t,
                                const BlockVector<InVectorType>& sol) override;

    virtual void postProcess(const double& t,
                             const BlockVector<InVectorType>& sol) override;

    virtual BlockMatrix<InMatrixType> getMass(const double& time,
                                      const BlockVector<InVectorType>& sol) override;

    virtual BlockMatrix<InMatrixType> getMassJacobian(const double& time,
                                      const BlockVector<InVectorType>& sol) override;

    virtual BlockVector<InVectorType> getRightHandSide(const double& time,
                                      const BlockVector<InVectorType>& sol) override;

    virtual BlockMatrix<InMatrixType> getJacobianRightHandSide(const double& time,
                                      const BlockVector<InVectorType>& sol) override;

    virtual BlockVector<InVectorType> getLifting(const double& time) const override;

    virtual BlockVector<InVectorType> getZeroVector() const override;

    virtual void apply0DirichletBCsMatrix(BlockMatrix<InMatrixType>& matrix, double diagCoeff) const override;

    virtual void apply0DirichletBCs(BlockVector<InVectorType>& vector) const override;

    virtual void setExporter() override;

    virtual void applyDirichletBCs(const double& time, BlockVector<InVectorType>& vector) const override;

    virtual void checkStabTerm(const BlockVector<InVectorType>& sol) const;

    std::map<unsigned int, std::string> getIDMeshTypeMap() const;

    inline SHP(aAssembler<VInner COMMA MInner>) block(const unsigned int& index) {return M_primalAssemblers[index];}

    std::map<unsigned int, SHP(aAssembler<VInner COMMA MInner>)> getAssemblersMap() const {return M_primalAssemblers;}

    std::vector<SHP(InterfaceAssembler<VInner COMMA MInner>)> getDualAssemblers() const {return M_dualAssemblers;}

    BlockVector<BlockVector<VectorEp>> convertFunctionRBtoFEM(BlockVector<BlockVector<DenseVector>> rbFunction, EPETRACOMM comm) const;

    virtual void setExtrapolatedSolution(const BlockVector<InVectorType>& exSol) override;

    virtual BlockVector<InVectorType> getNonLinearTerm() override;

    std::map<unsigned int,std::vector<double>> getRandomizibleParametersVectors();

    virtual void applyPiola(BlockVector<FEVECTOR> solution, bool inverse) override {};

    void applyGlobalPiola(BlockVector<BlockVector<FEVECTOR>> solution, bool inverse);

    virtual void initializeFEspaces() override;

protected:
    GetPot                                                        M_datafile;
    TreeStructure                                                 M_tree;
    std::map<unsigned int, SHP(aAssembler<VInner COMMA MInner>)>  M_primalAssemblers;
    std::vector<SHP(InterfaceAssembler<VInner COMMA MInner>)>     M_dualAssemblers;
    unsigned int                                                  M_numberBlocks;
    SHP(MDEIMManager)                                             M_mdeimManager;
    SHP(RBBasesManager)                                           M_basesManager;
};

}

#include "BlockAssembler_imp.hpp"

#endif // BLOCKASSEMBLER_HPP
