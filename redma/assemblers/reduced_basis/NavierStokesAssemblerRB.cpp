#include "NavierStokesAssemblerRB.hpp"

namespace RedMA
{


NavierStokesAssemblerRB::
NavierStokesAssemblerRB(const DataContainer& data, SHP(TreeNode) treeNode) :
  StokesAssemblerRB(data,treeNode),
  NavierStokesModel(data,treeNode)
{
    M_name = "NavierStokesAssemblerRB";
}

void
NavierStokesAssemblerRB::
addConvectiveMatrixRightHandSide(SHP(aVector) sol, SHP(aMatrix) mat)
{
    // note: the jacobian is not consistent but we do this for efficiency
}

void
NavierStokesAssemblerRB::
addConvectiveTermJacobianRightHandSide(SHP(aVector) sol, SHP(aVector) lifting,
                                       SHP(aMatrix) mat)
{
    // note: the jacobian is not consistent but we do this for efficiency
}

SHP(aVector)
NavierStokesAssemblerRB::
getRightHandSide(const double& time, const SHP(aVector)& sol)
{
    using namespace LifeV;
    using namespace ExpressionAssembly;

    Chrono chrono;
    chrono.start();

    std::string msg = "[NavierStokesAssembler] computing rhs ...";
    printlog(YELLOW, msg, this->M_data.getVerbose());

    SHP(aVector) retVec = StokesAssemblerRB::getRightHandSide(time,sol);

    bool approximatenonlinearterm = M_data("rb/online/approximatenonlinearterm", 1);

    if (!approximatenonlinearterm)
    {
        SHP(VECTOREPETRA)  nonLinearTerm(new VECTOREPETRA(M_velocityFESpace->map()));
        SHP(VECTOREPETRA)  velocityReconstructed;

        velocityReconstructed = M_bases->reconstructFEFunction(sol->block(0), 0, StokesModel::M_treeNode->M_ID);

        integrate(elements(M_velocityFESpaceETA->mesh()),
                   M_velocityFESpace->qr(),
                   M_velocityFESpaceETA,
                   value(this->M_density) *
                   dot(value(M_velocityFESpaceETA , *velocityReconstructed) *
                       grad(M_velocityFESpaceETA , *velocityReconstructed),
                   phi_i)
                 ) >> nonLinearTerm;

        nonLinearTerm->globalAssemble();

        SHP(DistributedVector) nonLinearTermWrap(new DistributedVector());
        nonLinearTermWrap->setVector(nonLinearTerm);

        SHP(BlockVector) nonLinearTermVec(new BlockVector(2));
        nonLinearTermVec->setBlock(0,nonLinearTermWrap);

        M_bcManager->apply0DirichletBCs(*nonLinearTermVec,
                                        getFESpaceBCs(),
                                        getComponentBCs());

        M_nonLinearTerm = M_bases->leftProject(nonLinearTermVec, ID());
    }
    else
    {
        unsigned int nterms = M_nonLinearTermsDecomposition.size();
        M_nonLinearTerm->multiplyByScalar(0);
        for (unsigned int i = 0; i < nterms; i++)
        {
            for (unsigned int j = 0; j < nterms; j++)
            {
                SHP(BlockVector) currVec(new BlockVector(0));
                currVec->hardCopy(M_nonLinearTermsDecomposition[i][j]);
                currVec->multiplyByScalar(sol->block(0)->operator()(i) *
                                          sol->block(0)->operator()(j));
                M_nonLinearTerm->add(currVec);
            }
        }
    }

    M_nonLinearTerm->multiplyByScalar(-1);
    // this->addNeumannBCs(retVec, time, sol);
    // retVec->add(M_nonLinearTerm);
    retVec->block(0)->add(M_nonLinearTerm->block(0));
    M_nonLinearTerm->multiplyByScalar(-1);

    msg = "done, in ";
    msg += std::to_string(chrono.diff());
    msg += " seconds\n";
    printlog(YELLOW, msg, this->M_data.getVerbose());

    return retVec;
}

SHP(aMatrix)
NavierStokesAssemblerRB::
getJacobianRightHandSide(const double& time,
                         const SHP(aVector)& sol)
{
    return StokesAssemblerRB::getJacobianRightHandSide(time,sol);
}

void
NavierStokesAssemblerRB::
RBsetup()
{
    using namespace LifeV;
    using namespace ExpressionAssembly;

    StokesAssemblerRB::RBsetup();

    if (M_data("rb/online/approximatenonlinearterm",1))
    {
        printlog(YELLOW, "[NavierStokesAssembler] precomputing non linear terms \t", M_data.getVerbose());
        Chrono chrono;
        chrono.start();

        auto velocityBasis = M_bases->getEnrichedBasis(0, StokesModel::M_treeNode->M_ID);

        int nterms = M_data("rb/online/numbernonlinearterms", 20);

        if (nterms == -1)
            nterms = velocityBasis.size();

        SHP(MATRIXEPETRA) nonLinearMatrix(new MATRIXEPETRA(M_velocityFESpace->map()));
        SHP(VECTOREPETRA) nonLinearTerm(new VECTOREPETRA(M_velocityFESpace->map()));

        M_nonLinearTermsDecomposition.resize(nterms);
        for (unsigned int i = 0; i < nterms; i++)
        {
            M_nonLinearTermsDecomposition[i].resize(nterms);
            nonLinearMatrix->zero();

            integrate(elements(M_velocityFESpaceETA->mesh()),
                       M_velocityFESpace->qr(),
                       M_velocityFESpaceETA,
                       M_velocityFESpaceETA,
                       value(this->M_density) *
                       dot(value(M_velocityFESpaceETA , *velocityBasis[i]) * grad(phi_j),
                       phi_i)
                     ) >> nonLinearMatrix;

            nonLinearMatrix->globalAssemble();

            for (unsigned int j = 0; j < nterms; j++)
            {
                *nonLinearTerm = (*nonLinearMatrix) * (*velocityBasis[j]);

                SHP(DistributedVector) nonLinearTermWrap(new DistributedVector());
                nonLinearTermWrap->setVector(nonLinearTerm);

                SHP(BlockVector) nonLinearTermVec(new BlockVector(2));
                nonLinearTermVec->setBlock(0,nonLinearTermWrap);

                M_bcManager->apply0DirichletBCs(*nonLinearTermVec,
                                                getFESpaceBCs(),
                                                getComponentBCs());
                M_nonLinearTermsDecomposition[i][j] = M_bases->leftProject(nonLinearTermVec, ID());
            }
        }

        std::string msg = "done, in ";
        msg += std::to_string(chrono.diff());
        msg += " seconds\n";
        printlog(YELLOW, msg, this->M_data.getVerbose());
    }
}

}
