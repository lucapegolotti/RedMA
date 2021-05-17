#include "BlockAssembler.hpp"

namespace RedMA
{
//
BlockAssembler::
BlockAssembler(const DataContainer& data,
               const TreeStructure& tree,
               shp<DefaultAssemblers> defAssemblers) :
  aAssembler(data),
  M_tree(tree)
{
    this->M_defaultAssemblers = defAssemblers;
    this->setup();
}

void
BlockAssembler::
checkStabTerm(const shp<aVector>& sol) const
{
    // for (auto as: M_dualAssemblers)
    //     as->checkStabilizationTerm(sol, M_primalAssemblers.size());
}

shp<aVector>
BlockAssembler::
getLifting(const double& time) const
{
    shp<BlockVector> retVec(new BlockVector(M_numberBlocks));

    for (auto as : M_primalAssemblers)
        retVec->setBlock(as.first,as.second->getLifting(time));

    return retVec;
}

void
BlockAssembler::
setDefaultAssemblers(shp<DefaultAssemblers> defAssemblers)
{
    this->M_defaultAssemblers = defAssemblers;
    for (auto as : M_primalAssemblers)
        as.second->setDefaultAssemblers(this->M_defaultAssemblers);
}

void
BlockAssembler::
applyPiola(shp<aVector> solution,
           bool inverse)
{
    unsigned int count = 0;
    for (auto as : M_primalAssemblers)
    {
        as.second->applyPiola(
            convert<BlockVector>(convert<BlockVector>(solution)->block(count)),
            inverse);
        count = count + 1;
    }
}

void
BlockAssembler::
applyDirichletBCsMatrix(shp<aMatrix> matrix,
                        double diagCoeff) const
{
    for (auto as : M_primalAssemblers)
        as.second->applyDirichletBCsMatrix(
            convert<BlockMatrix>(matrix)->block(as.first, as.first),
            diagCoeff);
}

void
BlockAssembler::
apply0DirichletBCs(shp<aVector> initialGuess) const
{
    for (auto as : M_primalAssemblers)
    {
        as.second->apply0DirichletBCs(
            convert<BlockVector>(initialGuess)->block(as.first));
    }
}

void
BlockAssembler::
applyDirichletBCs(const double& time,
                  shp<aVector> initialGuess) const
{
    for (auto as : M_primalAssemblers)
        as.second->applyDirichletBCs(time,
            convert<BlockVector>(initialGuess)->block(as.first));
}

shp<aVector>
BlockAssembler::
getZeroVector() const
{
    shp<BlockVector> retVec(new BlockVector(M_numberBlocks));

    for (auto as : M_primalAssemblers)
        retVec->setBlock(as.first,as.second->getZeroVector());

    unsigned int count = M_primalAssemblers.size();

    for (auto as : M_dualAssemblers)
    {
        auto newvec = as->getZeroVector();
        retVec->setBlock(count, newvec);
        count++;
    }
    // retVec->close();
    return retVec;
}

void
BlockAssembler::
setExporter()
{
    for (auto as: M_primalAssemblers)
        as.second->setExporter();
}

void
BlockAssembler::
exportSolution(const double& t, const shp<aVector>& sol)
{
    for (auto as : M_primalAssemblers)
        as.second->exportSolution(t,
            convert<BlockVector>(sol)->block(as.first));
}

std::map<unsigned int,std::vector<double>>
BlockAssembler::
getRandomizibleParametersVectors()
{
    std::map<unsigned int,std::vector<double>> retMap;

    for (auto as : M_primalAssemblers)
        retMap[as.first] = as.second->getTreeNode()->M_block->
                getGeometricParametersHandler().getRandomizibleParametersValueAsVector();

    return retMap;
}

void
BlockAssembler::
setExtrapolatedSolution(const shp<aVector>& exSol)
{
    for (auto as : M_primalAssemblers)
        as.second->setExtrapolatedSolution(
            convert<BlockVector>(exSol)->block(as.first));
}

void
BlockAssembler::
postProcess(const double& t, const shp<aVector>& sol)
{
    for (auto as : M_primalAssemblers)
        as.second->postProcess(t, convert<BlockVector>(sol)->block(as.first));

    if (this->M_data("coupling/check_stabterm", false))
        checkStabTerm(sol);
}

shp<aMatrix>
BlockAssembler::
getMass(const double& time,
        const shp<aVector>& sol)
{
    shp<BlockMatrix> mass(new BlockMatrix(M_numberBlocks, M_numberBlocks));

    for (auto as : M_primalAssemblers)
    {
        unsigned int ind = as.first;
        mass->setBlock(ind, ind, as.second->getMass(time, convert<BlockVector>(sol)->block(ind)));
    }

    return mass;
}

shp<aMatrix>
BlockAssembler::
getPressureMass(const double& time,
                const shp<aVector>& sol)
{
    shp<BlockMatrix> mass_press(new BlockMatrix(M_numberBlocks, M_numberBlocks));

    for (auto as : M_primalAssemblers)
    {
        unsigned int ind = as.first;
        mass_press->setBlock(ind, ind, as.second->getPressureMass(time, convert<BlockVector>(sol)->block(ind)));
    }

    return mass_press;
}

shp<aMatrix>
BlockAssembler::
getMassJacobian(const double& time,
                const shp<aVector>& sol)
{
    shp<BlockMatrix> massJacobian(new BlockMatrix(M_numberBlocks, M_numberBlocks));

    for (auto as : M_primalAssemblers)
    {
        unsigned int ind = as.first;
        massJacobian->setBlock(ind,ind,as.second->getMassJacobian(time, convert<BlockVector>(sol)->block(ind)));
    }

    return massJacobian;
}

shp<aVector>
BlockAssembler::
getRightHandSide(const double& time,
                 const shp<aVector>& sol)
{
    shp<BlockVector> rhs(new BlockVector(M_numberBlocks));

    for (auto as: M_primalAssemblers)
    {
        unsigned int ind = as.first;
        rhs->setBlock(ind, as.second->getRightHandSide(time, convert<BlockVector>(sol)->block(ind)));
    }

    // add interface contributions
    for (auto as: M_dualAssemblers)
    {
        as->addContributionRhs(time, rhs, convert<BlockVector>(sol),
            M_primalAssemblers.size());
    }

    return rhs;
}

shp<aVector>
BlockAssembler::
convertFunctionRBtoFEM(shp<aVector> rbFunction,
                       EPETRACOMM comm) const
{
    shp<BlockVector> retVec(new BlockVector(rbFunction->nRows()));

    for (auto as : M_primalAssemblers)
    {
        unsigned int ind = as.first;
        retVec->setBlock(ind,
            spcast<aAssemblerRB>(as.second)->convertFunctionRBtoFEM(convert<BlockVector>(rbFunction)->block(ind)));
    }

    if (rbFunction->nRows() > M_primalAssemblers.size())
    {
        for (auto as : M_dualAssemblers)
        {
            unsigned int indInterface = as->getInterface().M_ID + M_primalAssemblers.size();
            spcast<BlockVector>(retVec->block(indInterface))->resize(1);

            convert<BlockVector>(retVec->block(indInterface))->setBlock(0,
            DistributedVector::convertDenseVector(
            convert<DenseVector>(
                convert<BlockVector>(
                    convert<BlockVector>(rbFunction)
                    ->block(indInterface))->block(0)),comm));
        }
    }
    return retVec;
}

shp<aVector>
BlockAssembler::
getNonLinearTerm()
{
    shp<BlockVector> retVec(new BlockVector(M_primalAssemblers.size()));

    for (auto as : M_primalAssemblers)
        retVec->setBlock(as.first,as.second->getNonLinearTerm());

    return retVec;
}


shp<aMatrix>
BlockAssembler::
getJacobianRightHandSide(const double& time,
                         const shp<aVector>& sol)
{
    shp<BlockMatrix> jac(new BlockMatrix(M_numberBlocks, M_numberBlocks));

    for (auto as: M_primalAssemblers)
    {
        unsigned int ind = as.first;
        jac->setBlock(ind,ind,as.second->getJacobianRightHandSide(time,
                              convert<BlockVector>(sol)->block(ind)));
    }

    for (auto as: M_dualAssemblers)
        as->addContributionJacobianRhs(time, jac,
            spcast<BlockVector>(sol),
             M_primalAssemblers.size());

    return jac;
}

std::map<unsigned int, std::string>
BlockAssembler::
getIDMeshTypeMap() const
{
    std::map<unsigned int, std::string> retMap;
    for (auto as: M_primalAssemblers)
    {
        std::string meshname = as.second->getTreeNode()->M_block->getMeshName();
        unsigned int dashpos = meshname.find("/");
        unsigned int formatpos = meshname.find(".mesh");
        std::string actualmeshname = meshname.substr(dashpos + 1, formatpos - dashpos - 1);
        retMap[as.first] = actualmeshname;
    }

    return retMap;
}

void
BlockAssembler::
setup()
{
    typedef std::map<unsigned int, shp<TreeNode>>         NodesMap;
    typedef aAssembler                                    InnerAssembler;
    typedef std::vector<shp<TreeNode>>                    NodesVector;

    printlog(GREEN, "[BlockAssembler] initializing block assembler ... \n", this->M_data.getVerbose());

    NodesMap nodesMap = M_tree.getNodesMap();
    // allocate assemblers
    for (auto it = nodesMap.begin(); it != nodesMap.end(); it++)
    {
        shp<InnerAssembler> newAssembler;
        newAssembler = AssemblerFactory(this->M_data, it->second);
        newAssembler->setDefaultAssemblers(M_defaultAssemblers);
        M_comm = newAssembler->getComm();
        M_primalAssemblers[it->second->M_ID] = newAssembler;
    }

    if (!(this->arePrimalAssemblersFE())) {
        M_basesManager.reset(new RBBasesManager(M_data, M_comm, getIDMeshTypeMap()));
        M_basesManager->load();
    }

    for (auto& primalas : M_primalAssemblers)
    {
        primalas.second->setup();
        // we set the RBBases later because we need the finite element space
        auto block = primalas.second->getTreeNode()->M_block;
        if (std::strcmp(block->getDiscretizationMethod().c_str(), "rb") == 0)
            primalas.second->setRBBases(M_basesManager);
    }

    // we need to load the bases here because now we have the finite element space
    if(!(this->arePrimalAssemblersFE()))
        M_basesManager->loadBases();

    for (auto& primalas : M_primalAssemblers)
    {
        // restrict RB matrices based on desired pod tolerance (if needed)
        auto block = primalas.second->getTreeNode()->M_block;
        if (std::strcmp(block->getDiscretizationMethod().c_str(), "rb") == 0)
            primalas.second->RBsetup();
    }

    // allocate interface assemblers
    unsigned int interfaceID = 0;
    for (NodesMap::iterator it = nodesMap.begin(); it != nodesMap.end(); it++)
    {
        NodesVector children = it->second->M_children;

        unsigned int myID = it->second->M_ID;

        shp<InnerAssembler> fatherAssembler = M_primalAssemblers[myID];

        unsigned int countChildren = 0;
        for (NodesVector::iterator itVector = children.begin();
             itVector != children.end(); itVector++)
        {
            if (*itVector)
            {
                unsigned int otherID = (*itVector)->M_ID;
                shp<InnerAssembler> childAssembler = M_primalAssemblers[otherID];

                if (spcast<StokesAssemblerFE>(fatherAssembler)->hasNoSlipBCs() !=
                    spcast<StokesAssemblerFE>(childAssembler)->hasNoSlipBCs())
                    throw new Exception("Father and Child assemblers MUST have the same "
                                        "type of BCs at the vessel wall!");

                Interface newInterface(fatherAssembler, myID,
                                       childAssembler, otherID,
                                       interfaceID);
                newInterface.M_indexOutlet = countChildren;
                shp<InterfaceAssembler> inAssembler;
                inAssembler.reset(new InterfaceAssembler(this->M_data, newInterface,
                                                           spcast<StokesAssemblerFE>(fatherAssembler)->hasNoSlipBCs()));
                M_dualAssemblers.push_back(inAssembler);
                interfaceID++;
            }
            countChildren++;
        }
    }

    if ((!std::strcmp(this->M_data("bc_conditions/inlet_bc_type", "dirichlet").c_str(), "dirichlet")) &&
        (!std::strcmp(this->M_data("bc_conditions/inletdirichlet", "weak").c_str(), "weak")))
    {
        shp<InnerAssembler> inletAssembler = M_primalAssemblers[0];

        // we set the inlet to child such that we are consistent with the normal orientation
        // with respect to flow direction
        Interface newInterface(nullptr, -1, inletAssembler, 0,
                                               interfaceID);

        shp<InterfaceAssembler> inletInAssembler;
        bool inletPrimalIsFE = !(std::strcmp(inletAssembler->getTreeNode()->M_block->getDiscretizationMethod().c_str(), "fem"));
        bool hasNoSlipBCs = (inletPrimalIsFE) ? ((spcast<StokesAssemblerFE>(inletAssembler))->hasNoSlipBCs()) :
                            (spcast<StokesAssemblerFE>(spcast<StokesAssemblerRB>(inletAssembler)->getFEAssembler())->hasNoSlipBCs());
        inletInAssembler.reset(new InletInflowAssembler(this->M_data, newInterface, hasNoSlipBCs));

        M_dualAssemblers.push_back(inletInAssembler);
        interfaceID++;
    }

    M_numberBlocks = M_primalAssemblers.size() + M_dualAssemblers.size();

    printlog(GREEN, "done\n", this->M_data.getVerbose());
}

bool
BlockAssembler::
arePrimalAssemblersFE()
{
    for (auto& primalas : M_primalAssemblers)
    {
        auto block = primalas.second->getTreeNode()->M_block;
        if (std::strcmp(block->getDiscretizationMethod().c_str(), "fem") != 0)
            return false;
    }
    return true;
}

bool
BlockAssembler::
arePrimalAssemblersRB()
{
    for (auto& primalas : M_primalAssemblers)
    {
        auto block = primalas.second->getTreeNode()->M_block;
        if (std::strcmp(block->getDiscretizationMethod().c_str(), "rb") != 0)
            return false;
    }
    return true;
}

}
