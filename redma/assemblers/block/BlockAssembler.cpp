#include "BlockAssembler.hpp"

namespace RedMA
{

BlockAssembler::
BlockAssembler(const DataContainer& data, const TreeStructure& tree,
               SHP(DefaultAssemblers) defAssemblers) :
  aAssembler(data),
  M_tree(tree)
{
    // this->M_defaultAssemblers = defAssemblers;
    // setup();
}

void
BlockAssembler::
checkStabTerm(const BlockVector& sol) const
{
    // for (auto as: M_dualAssemblers)
    //     as->checkStabilizationTerm(sol, M_primalAssemblers.size());
}

void
BlockAssembler::
setExporter()
{
    // for (auto as: M_primalAssemblers)
    //     as.second->setExporter();
}

void
BlockAssembler::
initializeFEspaces()
{
    // for (auto as: M_primalAssemblers)
    //     as.second->initializeFEspaces();
}

BlockVector
BlockAssembler::
getLifting(const double& time) const
{
    // BlockVector retVec(M_numberBlocks);
    //
    // for (auto as : M_primalAssemblers)
    //     retVec.block(as.first) = as.second->getLifting(time);
    //
    // return retVec;
}

void
BlockAssembler::
setDefaultAssemblers(SHP(DefaultAssemblers) defAssemblers)
{
    // this->M_defaultAssemblers = defAssemblers;
    // for (auto as : M_primalAssemblers)
    //     as.second->setDefaultAssemblers(this->M_defaultAssemblers);
}

void
BlockAssembler::
applyGlobalPiola(BlockVector solution, bool inverse)
{
    // unsigned int count = 0;
    // for (auto as : M_primalAssemblers)
    // {
    //     as.second->applyPiola(solution.block(count), inverse);
    //     count = count + 1;
    // }
}

void
BlockAssembler::
apply0DirichletBCsMatrix(BlockMatrix& matrix, double diagCoeff) const
{
    // for (auto as : M_primalAssemblers)
    //     as.second->apply0DirichletBCsMatrix(matrix.block(as.first, as.first), diagCoeff);
}

void
BlockAssembler::
apply0DirichletBCs(BlockVector& initialGuess) const
{
    // for (auto as : M_primalAssemblers)
    //     as.second->apply0DirichletBCs(initialGuess.block(as.first));
}

void
BlockAssembler::
applyDirichletBCs(const double& time, BlockVector& initialGuess) const
{
    // for (auto as : M_primalAssemblers)
    //     as.second->applyDirichletBCs(time, initialGuess.block(as.first));
}

BlockVector
BlockAssembler::
getZeroVector() const
{
    // BlockVector retVec;
    // retVec.resize(M_numberBlocks);
    //
    // for (auto as : M_primalAssemblers)
    //     retVec.block(as.first).softCopy(as.second->getZeroVector());
    //
    // unsigned int count = M_primalAssemblers.size();
    // for (auto as : M_dualAssemblers)
    // {
    //     retVec.block(count).softCopy(as->getZeroVector());
    //     count++;
    // }
    //
    // return retVec;
}


void
BlockAssembler::
exportSolution(const double& t, const BlockVector& sol)
{
    // for (auto as : M_primalAssemblers)
    //     as.second->exportSolution(t, sol.block(as.first));
}

std::map<unsigned int,std::vector<double>>
BlockAssembler::
getRandomizibleParametersVectors()
{
    // std::map<unsigned int,std::vector<double>> retMap;
    //
    // for (auto as : M_primalAssemblers)
    //     retMap[as.first] = as.second->getTreeNode()->M_block->
    //             getGeometricParametersHandler().getRandomizibleParametersValueAsVector();
    //
    // return retMap;
}

void
BlockAssembler::
setExtrapolatedSolution(const BlockVector& exSol)
{
    // for (auto as : M_primalAssemblers)
    //     as.second->setExtrapolatedSolution(exSol.block(as.first));
}

void
BlockAssembler::
postProcess(const double& t, const BlockVector& sol)
{
    // for (auto as : M_primalAssemblers)
    //     as.second->postProcess(t, sol.block(as.first));
    //
    // if (this->M_data("coupling/check_stabterm", false))
    //     checkStabTerm(sol);

}

BlockMatrix
BlockAssembler::
getMass(const double& time, const BlockVector& sol)
{
    // BlockMatrix<InMatrixType> mass;
    // mass.resize(M_numberBlocks, M_numberBlocks);
    //
    // for (auto as : M_primalAssemblers)
    // {
    //     unsigned int ind = as.first;
    //     mass.block(ind, ind).softCopy(as.second->getMass(time, sol.block(ind)));
    // }
    //
    // return mass;
}

BlockMatrix
BlockAssembler::
getMassJacobian(const double& time, const BlockVector& sol)
{
    // BlockMatrix massJacobian;
    // massJacobian.resize(M_numberBlocks, M_numberBlocks);
    //
    // for (auto as : M_primalAssemblers)
    // {
    //     unsigned int ind = as.first;
    //     massJacobian.block(ind, ind).softCopy(as.second->getMassJacobian(time, sol.block(ind)));
    // }
    //
    // return massJacobian;
}

BlockVector
BlockAssembler::
getRightHandSide(const double& time, const BlockVector& sol)
{
    // BlockVector rhs;
    // rhs.resize(M_numberBlocks);
    //
    // for (auto as: M_primalAssemblers)
    // {
    //     unsigned int ind = as.first;
    //     rhs.block(ind).softCopy(as.second->getRightHandSide(time, sol.block(ind)));
    // }
    //
    // // add interface contributions
    // for (auto as: M_dualAssemblers)
    //     as->addContributionRhs(time, rhs, sol, M_primalAssemblers.size());
    //
    // return rhs;
}

BlockVector
BlockAssembler::
convertFunctionRBtoFEM(BlockVector rbFunction,
                       EPETRACOMM comm) const
{
    // BlockVector retVec(rbFunction.nRows());
    //
    // for (auto as : M_primalAssemblers)
    // {
    //     unsigned int ind = as.first;
    //     retVec.block(ind).softCopy(as.second->convertFunctionRBtoFEM(rbFunction.block(ind)));
    // }
    //
    // if (rbFunction.nRows() > M_primalAssemblers.size())
    // {
    //     for (auto as : M_dualAssemblers)
    //     {
    //         unsigned int indInterface = as->getInterface().M_ID + M_primalAssemblers.size();
    //         retVec.block(indInterface).resize(1);
    //         retVec.block(indInterface).block(0).softCopy(
    //                 VectorEp::convertDenseVector(rbFunction.block(indInterface).block(0),
    //                 comm));
    //     }
    // }
    // return retVec;
}

BlockVector
BlockAssembler::
getNonLinearTerm()
{
    // BlockVector retVec(M_primalAssemblers.size());
    //
    // for (auto as : M_primalAssemblers)
    //     retVec.block(as.first) = as.second->getNonLinearTerm();
    //
    // return retVec;
}


BlockMatrix
BlockAssembler::
getJacobianRightHandSide(const double& time, const BlockVector& sol)
{
    // BlockMatrix jac;
    // jac.resize(M_numberBlocks, M_numberBlocks);
    //
    // for (auto as: M_primalAssemblers)
    // {
    //     unsigned int ind = as.first;
    //     jac.block(ind,ind).softCopy(as.second->getJacobianRightHandSide(time,
    //                                 sol.block(ind)));
    // }
    //
    // for (auto as: M_dualAssemblers)
    //     as->addContributionJacobianRhs(time, jac, sol, M_primalAssemblers.size());
    //
    // return jac;
}

std::map<unsigned int, std::string>
BlockAssembler::
getIDMeshTypeMap() const
{
    // std::map<unsigned int, std::string> retMap;
    // for (auto as: M_primalAssemblers)
    // {
    //     std::string meshname = as.second->getTreeNode()->M_block->getMeshName();
    //     unsigned int dashpos = meshname.find("/");
    //     unsigned int formatpos = meshname.find(".mesh");
    //     std::string actualmeshname = meshname.substr(dashpos + 1, formatpos - dashpos - 1);
    //     retMap[as.first] = actualmeshname;
    // }
    //
    // return retMap;
}

}
