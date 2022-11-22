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
    /*for (auto as: M_dualAssemblers)
        as->checkStabilizationTerm(spcast<BlockVector>(sol),
                        M_primalAssemblers.size());*/
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

shp<aVector>
BlockAssembler::
getDisplacement() const
{
    shp<BlockVector> retVec(new BlockVector(M_numberBlocks));

    for (auto as : M_primalAssemblers)
        retVec->setBlock(as.first,as.second->getDisplacement());

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

    return retVec;
}

std::map<unsigned int, std::vector<shp<BlockVector>>>
BlockAssembler::
importSolution(const std::string& filename) const
{
    if (!fs::exists(filename))
        throw new Exception("Importing error. Invalid path provided!");

    printlog(GREEN, "[BlockAssembler] importing solution ...\n", this->M_data.getVerbose());

    unsigned int cnt = 0;
    std::map<unsigned int, std::vector<shp<BlockVector>>> retMap;
    for (auto as : M_primalAssemblers)
    {
        std::string fname = filename + "Block" + std::to_string(cnt) + "/";
        std::vector<shp<BlockVector>> sol = as.second->importSolution(fname)[0];
        retMap[cnt] = sol;
        cnt += 1;
    }

    return retMap;
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

    /*if (this->M_data("coupling/check_stabterm", false))
        checkStabTerm(sol);*/
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
        M_basesManager->loadSingularValues();
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
                unsigned int ninlets = childAssembler->getTreeNode()->M_block->getInlets().size();

                if (spcast<StokesAssemblerFE>(fatherAssembler)->hasNoSlipBCs() !=
                    spcast<StokesAssemblerFE>(childAssembler)->hasNoSlipBCs())
                    throw new Exception("Father and Child assemblers MUST have the same "
                                        "type of BCs at the vessel wall!");

                for (unsigned int i = 0; i < ninlets; i++)
                {
                    if ((childAssembler->getTreeNode()->M_block->getInlet(i) ==
                        fatherAssembler->getTreeNode()->M_block->getOutlet(countChildren)) ||
                        (ninlets == 1))
                    {
                        Interface newInterface(fatherAssembler, myID,
                                               childAssembler, otherID,
                                               interfaceID);
                        newInterface.M_indexInlet = i;
                        newInterface.M_indexOutlet = countChildren;
                        shp<InterfaceAssembler> inAssembler;
                        inAssembler.reset(new InterfaceAssembler(this->M_data, newInterface,
                                                                 spcast<StokesAssemblerFE>(fatherAssembler)->hasNoSlipBCs()));
                        M_dualAssemblers.push_back(inAssembler);
                        interfaceID++;

                        i = ninlets; // once I find a matching inlet in the child, I can avoid investigating the others
                    }
                }
            }
            countChildren++;
        }
    }

    // allocate weak inlet Dirichlet BC assemblers
    if ((!std::strcmp(this->M_data("bc_conditions/inlet_bc_type", "dirichlet").c_str(), "dirichlet")) &&
        (!std::strcmp(this->M_data("bc_conditions/inletdirichlet", "weak").c_str(), "weak")))
    {
        shp<InnerAssembler> inletAssembler = M_primalAssemblers[0];

        std::vector<GeometricFace> inlets = inletAssembler->getTreeNode()->M_block->getInlets();
        unsigned int ninlets = inlets.size();

        for (unsigned int i = 0; i < ninlets; i++)
        {
            // we set the inlet to child such that we are consistent with the normal orientation
            // with respect to flow direction
            Interface newInterface(nullptr, -1, inletAssembler, 0,
                                                   interfaceID);
            newInterface.M_indexInlet = i;
            newInterface.M_interfaceFlag = inlets[i].M_flag;

            shp<InterfaceAssembler> inletInAssembler;
            bool inletPrimalIsFE = !(std::strcmp(inletAssembler->getTreeNode()->M_block->getDiscretizationMethod().c_str(), "fem"));
            bool hasNoSlipBCs = (inletPrimalIsFE) ? ((spcast<StokesAssemblerFE>(inletAssembler))->hasNoSlipBCs()) :
                            (spcast<StokesAssemblerFE>(spcast<StokesAssemblerRB>(inletAssembler)->getFEAssembler())->hasNoSlipBCs());
            inletInAssembler.reset(new InletInflowAssembler(this->M_data, newInterface, hasNoSlipBCs));

            M_dualAssemblers.push_back(inletInAssembler);
            interfaceID++;
        }
    }

    // allocate weak outlet Dirichlet BC assemblers
    unsigned int numConditions = this->M_data("bc_conditions/numoutletbcs", 0);
    for (unsigned int outletIndex = 0; outletIndex < numConditions; outletIndex++)
    {
        std::string dataEntry = "bc_conditions/outlet" + std::to_string(outletIndex);
        unsigned int blockindex = this->M_data(dataEntry + "/blockindex", 0);
        std::string BCtype = this->M_data(dataEntry + "/type", "windkessel");

        if (!std::strcmp(BCtype.c_str(), "dirichlet"))
        {
            shp<InnerAssembler> outletAssembler = M_primalAssemblers[blockindex];
            auto outlets = outletAssembler->getTreeNode()->M_block->getOutlets();
            unsigned int boundaryflag = this->M_data(dataEntry + "/boundaryflag", 2);

            for (GeometricFace outlet : outlets)
            {
                unsigned int cnt = 0;
                if (outlet.M_flag == boundaryflag)
                {
                    Interface newInterface(outletAssembler, blockindex,
                                           nullptr, -1, interfaceID);
                    newInterface.M_indexOutlet = cnt;
                    newInterface.M_interfaceFlag = outlet.M_flag;

                    shp<OutletOutflowAssembler> outletOutAssembler;
                    bool outletPrimalIsFE = !(std::strcmp(outletAssembler->getTreeNode()->M_block->getDiscretizationMethod().c_str(), "fem"));
                    bool hasNoSlipBCs = (outletPrimalIsFE) ? ((spcast<StokesAssemblerFE>(outletAssembler))->hasNoSlipBCs()) :
                            (spcast<StokesAssemblerFE>(spcast<StokesAssemblerRB>(outletAssembler)->getFEAssembler())->hasNoSlipBCs());
                    outletOutAssembler.reset(new OutletOutflowAssembler(this->M_data, newInterface, hasNoSlipBCs));
                    outletOutAssembler->setGlobalOutletIndex(outletIndex);

                    M_dualAssemblers.push_back(outletOutAssembler);
                    interfaceID++;
                }
                ++cnt;
            }
        }
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
