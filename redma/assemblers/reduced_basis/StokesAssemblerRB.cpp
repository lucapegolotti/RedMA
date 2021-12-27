#include "StokesAssemblerRB.hpp"

namespace RedMA
{

StokesAssemblerRB::
StokesAssemblerRB(const DataContainer& data, shp<TreeNode> treeNode) :
  aAssemblerRB(data, treeNode)
{
    M_FEAssembler.reset(new StokesAssemblerFE(data, treeNode));
    M_name = "StokesAssemblerRB";
    M_nComponents = 2;
}

void
StokesAssemblerRB::
setup()
{
    std::string msg = "[";
    msg += this->M_name;
    msg += "] initializing internal FE assembler...";
    printlog(YELLOW, msg, this->M_data.getVerbose());
    M_FEAssembler->setup();
}

shp<aMatrix>
StokesAssemblerRB::
getMass(const double& time,
        const shp<aVector>& sol)
{
    return M_reducedMass;
}

shp<aMatrix>
StokesAssemblerRB::
getPressureMass(const double& time,
                const shp<aVector>& sol)
{
    return M_reducedMassPressure;
}

shp<aMatrix>
StokesAssemblerRB::
getMassJacobian(const double& time,
                const shp<aVector>& sol)
{
    shp<BlockMatrix> retMat(new BlockMatrix(this->M_nComponents,
                                              this->M_nComponents));
    return retMat;
}

shp<aVector>
StokesAssemblerRB::
getRightHandSide(const double& time,
                 const shp<aVector>& sol)
{
    shp<BlockMatrix> systemMatrix(new BlockMatrix(this->M_nComponents,
                                                  this->M_nComponents));
    systemMatrix->add(M_reducedStiffness);
    systemMatrix->add(M_reducedDivergence);
    systemMatrix->multiplyByScalar(-1.0);
    shp<aVector> retVec = systemMatrix->multiplyByVector(sol);
    return retVec;
}

shp<aMatrix>
StokesAssemblerRB::
getJacobianRightHandSide(const double& time,
                         const shp<aVector>& sol)
{
    shp<BlockMatrix> retMat(new BlockMatrix(this->M_nComponents,
                                            this->M_nComponents));

    retMat->add(M_reducedStiffness);
    retMat->add(M_reducedDivergence);
    retMat->multiplyByScalar(-1.0);

    return retMat;
}

std::vector<shp<aMatrix>>
StokesAssemblerRB::
getMatrices() const
{
    std::vector<shp<aMatrix>> retVec;

    retVec.push_back(M_reducedMass);
    retVec.push_back(M_reducedStiffness);
    retVec.push_back(M_reducedDivergence);

    return retVec;
}

shp<aMatrix>
StokesAssemblerRB::
assembleMatrix(const unsigned int& index)
{
    // throw new Exception("assembleMatrix not implemented for StokesAssemblerRB");
}

void
StokesAssemblerRB::
applyDirichletBCsMatrix(shp<aMatrix> matrix,
                        double diagCoeff) const
{
    // throw new Exception("Method not implemented for RB");
}

void
StokesAssemblerRB::
applyDirichletBCs(const double& time,
                  shp<aVector> vector) const
{
    // throw new Exception("Method not implemented for RB");
}

void
StokesAssemblerRB::
postProcess(const double& t, const shp<aVector>& sol)
{
    getBCManager()->postProcess();
}

void
StokesAssemblerRB::
apply0DirichletBCs(shp<aVector> vector) const
{
    // throw new Exception("Method not implemented for RB");
}

std::map<unsigned int, std::vector<shp<BlockVector>>>
StokesAssemblerRB::
importSolution(const std::string &filename) const
{
    return M_FEAssembler->importSolution(filename);
}

void
StokesAssemblerRB::
setExporter()
{
    this->M_FEAssembler->setExporter();
}

void
StokesAssemblerRB::
exportSolution(const double& t, 
               const shp<aVector>& sol)
{
    std::ofstream outfile("rbcoefs/block" + std::to_string(M_treeNode->M_ID) + "t" + std::to_string(t) + ".txt");
    std::string str2write = spcast<DenseVector>(sol->block(0))->getString(',') + "\n";
    outfile.write(str2write.c_str(), str2write.size());
    outfile.close();

    unsigned int id = M_treeNode->M_ID;

    shp<BlockVector> solBlock(new BlockVector(2));
    solBlock->setBlock(0,wrap(M_bases->reconstructFEFunction(sol->block(0), 0, id)));
    solBlock->setBlock(1,wrap(M_bases->reconstructFEFunction(sol->block(1), 1, id)));

    M_FEAssembler->exportSolution(t, solBlock);
}

shp<aVector>
StokesAssemblerRB::
getZeroVector() const
{
    shp<BlockVector> retVec(new BlockVector(M_nComponents));

    shp<DENSEVECTOR> uComp(new DENSEVECTOR(M_bases->getSizeEnrichedBasis(0)));
    uComp->Scale(0.0);
    shp<DENSEVECTOR> pComp(new DENSEVECTOR(M_bases->getSizeEnrichedBasis(1)));
    pComp->Scale(0.0);

    shp<DenseVector> uCompWrap(new DenseVector());
    uCompWrap->setVector(uComp);

    shp<DenseVector> pCompWrap(new DenseVector());
    pCompWrap->setVector(pComp);

    retVec->setBlock(0,uCompWrap);
    retVec->setBlock(1,pCompWrap);

    return retVec;
}

shp<aVector>
StokesAssemblerRB::
getLifting(const double& time) const
{
    return M_FEAssembler->getLifting(time);
}

void
StokesAssemblerRB::
RBsetup()
{
    if (M_bases == nullptr)
        throw new Exception("RB bases have not been set yet");

    // scale bases with piola transformation
    unsigned int indexField = 0;
    printlog(YELLOW, "[StokesAssemblerRB] applying Piola transformation\t", M_data.getVerbose());
    Chrono chrono;
    chrono.start();
    M_bases->scaleBasisWithPiola(0, M_treeNode->M_ID, [=](shp<VECTOREPETRA> vector)
    {
        shp<BlockVector> vectorWrap(new BlockVector(2));

        shp<DistributedVector> compVec(new DistributedVector());
        compVec->setVector(vector);

        vectorWrap->setBlock(0,compVec);

        applyPiola(vectorWrap, false);
    });
    std::string msg = "done, in ";
    msg += std::to_string(chrono.diff());
    msg += " seconds\n";
    printlog(YELLOW, msg, this->M_data.getVerbose());

    printlog(YELLOW, "[StokesAssemblerRB] assembling and projecting matrices\t", M_data.getVerbose());
    chrono.start();

    unsigned int id = M_treeNode->M_ID;

    M_reducedMass.reset(new BlockMatrix(2,2));
    M_reducedStiffness.reset(new BlockMatrix(2,2));
    M_reducedDivergence.reset(new BlockMatrix(2,2));
    M_reducedMassPressure.reset(new BlockMatrix(2,2));

    auto matrices = M_FEAssembler->getMatrices();
    M_reducedMass->setBlock(0,0,M_bases->matrixProject(matrices[0]->block(0,0), 0, 0, id));
    M_reducedStiffness->setBlock(0,0,M_bases->matrixProject(matrices[1]->block(0,0), 0, 0, id));
    M_reducedDivergence->setBlock(0,1,M_bases->matrixProject(matrices[2]->block(0,1), 0, 1, id));
    M_reducedDivergence->setBlock(1,0,M_bases->matrixProject(matrices[2]->block(1,0), 1, 0, id));

    auto mass_pressure = M_FEAssembler->getPressureMass(0.0, nullptr);
    M_reducedMassPressure->setBlock(1,1,M_bases->matrixProject(mass_pressure->block(1,1), 1, 1, id));

    msg = "done, in ";
    msg += std::to_string(chrono.diff());
    msg += " seconds\n";
    printlog(YELLOW, msg, this->M_data.getVerbose());
}

shp<RBBases>
StokesAssemblerRB::
getRBBases() const
{
    return M_bases;
}

void
StokesAssemblerRB::
setRBBases(shp<RBBasesManager> rbManager)
{
    // std::string mdeimdir = this->M_data("rb/online/mdeim/directory", "mdeims");
    std::string meshName = M_treeNode->M_block->getMeshName();
    unsigned int dashpos = meshName.find("/");
    unsigned int formatpos = meshName.find(".mesh");
    std::string actualName = meshName.substr(dashpos + 1,
                                             formatpos - dashpos - 1);

    // beware that at this point the rb bases have not been loaded yet
    M_bases = rbManager->getRBBases(actualName);
    M_bases->setFESpace(M_FEAssembler->getFEspace(0), 0);
    M_bases->setFESpace(M_FEAssembler->getFEspace(1), 1);
}

shp<aVector>
StokesAssemblerRB::
convertFunctionRBtoFEM(shp<aVector> rbSolution) const
{
    auto rbSolutionBlck = convert<BlockVector>(rbSolution);

    shp<BlockVector> retVec(new BlockVector(2));

    unsigned int id = M_treeNode->M_ID;

    if (rbSolutionBlck->block(0)->data())
    {
        shp<DenseVector> comp0 = convert<DenseVector>(rbSolutionBlck->block(0));
        shp<DistributedVector> rec0(new DistributedVector());
        rec0->setVector(M_bases->reconstructFEFunction(comp0, 0, id));
        retVec->setBlock(0,rec0);
    }

    if (rbSolutionBlck->block(1)->data())
    {
        shp<DenseVector> comp1 = convert<DenseVector>(rbSolutionBlck->block(1));
        shp<DistributedVector> rec1(new DistributedVector());
        rec1->setVector(M_bases->reconstructFEFunction(comp1, 1, id));
        retVec->setBlock(1,rec1);
    }

    return retVec;
}

void
StokesAssemblerRB::
applyPiola(shp<aVector> solution,
           bool inverse)
{
    M_FEAssembler->applyPiola(solution, inverse);
}

void
StokesAssemblerRB::
setDefaultAssemblers(shp<DefaultAssemblersLibrary> defAssemblers)
{
    M_defaultAssemblers = defAssemblers;
    M_FEAssembler->setDefaultAssemblers(defAssemblers);
}

shp<BlockVector>
StokesAssemblerRB::
reconstructFESolution(shp<BlockVector> sol) const
{
    shp<VECTOREPETRA>  velocityHandlerSol;
    shp<VECTOREPETRA>  pressureHandlerSol;
    velocityHandlerSol = M_bases->reconstructFEFunction(convert<BlockVector>(sol)->block(0),
                                                        0, ID());
    pressureHandlerSol = M_bases->reconstructFEFunction(convert<BlockVector>(sol)->block(1),
                                                        1, ID());

    shp<BlockVector> solFEM(new BlockVector(this->M_nComponents));
    solFEM->setBlock(0, wrap(velocityHandlerSol));
    solFEM->setBlock(1, wrap(pressureHandlerSol));

    return solFEM;
}

}
