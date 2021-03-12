#include "NavierStokesAssemblerRB.hpp"

namespace RedMA
{


NavierStokesAssemblerRB::
NavierStokesAssemblerRB(const DataContainer& data, shp<TreeNode> treeNode) :
  StokesAssemblerRB(data,treeNode)
{
    M_name = "NavierStokesAssemblerRB";
    M_exactJacobian = M_data("rb/online/exactjacobian",0);
}

void
NavierStokesAssemblerRB::
addConvectiveMatrixRightHandSide(shp<aVector> sol, shp<aMatrix> mat)
{
    // note: the jacobian is not consistent but we do this for efficiency
}

void
NavierStokesAssemblerRB::
addConvectiveTermJacobianRightHandSide(shp<aVector> sol, shp<aVector> lifting,
                                       shp<aMatrix> mat)
{
    Chrono chrono;
    chrono.start();

    std::string msg = "[NavierStokesAssemblerRB] computing Jacobian rhs ...";
    printlog(YELLOW, msg, this->M_data.getVerbose());

    using namespace LifeV;
    using namespace ExpressionAssembly;

    shp<VECTOREPETRA>  velocityHandler;
    velocityHandler = M_bases->reconstructFEFunction(convert<BlockVector>(sol)->block(0), 0, M_treeNode->M_ID);

    // shp<VECTOREPETRA>  liftingHandler;
    // liftingHandler = M_bases->reconstructFEFunction(lifting->block(0), 0, StokesModel::M_treeNode->M_ID);

    shp<MATRIXEPETRA>  convectiveMatrix(new MATRIXEPETRA(M_FEAssembler->getFEspace(0)->map()));
    shp<VECTOREPETRA>  velocityRepeated(new VECTOREPETRA(*velocityHandler, Repeated));
    // shp<VECTOREPETRA>  liftingRepeated(new VECTOREPETRA(*liftingHandler, Repeated));

    shp<ETFESPACE3> velocityFESpaceETA = M_FEAssembler->getVelocityETFEspace();
    double density = M_FEAssembler->getDensity();
    integrate(elements(velocityFESpaceETA->mesh()),
               M_FEAssembler->getFEspace(0)->qr(),
               velocityFESpaceETA,
               velocityFESpaceETA,
               value(density) *
               dot(
               (
               value(velocityFESpaceETA , *velocityRepeated) * grad(phi_j) +
               phi_j * grad(velocityFESpaceETA , *velocityRepeated)
               ),
               phi_i)
             ) >> convectiveMatrix;

    convectiveMatrix->globalAssemble();

    shp<SparseMatrix> matrixWrap(new SparseMatrix());
    matrixWrap->setMatrix(convectiveMatrix);

    unsigned int id = M_treeNode->M_ID;
    auto projectedMat = M_bases->matrixProject(matrixWrap, 0, 0, id);

    projectedMat->multiplyByScalar(-1);
    convert<BlockMatrix>(mat)->block(0,0)->add(projectedMat);

    msg = "done, in ";
    msg += std::to_string(chrono.diff());
    msg += " seconds\n";
    printlog(YELLOW, msg, this->M_data.getVerbose());
}

shp<aVector>
NavierStokesAssemblerRB::
getRightHandSide(const double& time, const shp<aVector>& sol)
{
    using namespace LifeV;
    using namespace ExpressionAssembly;

    Chrono chrono;
    chrono.start();

    auto solBlck = convert<BlockVector>(sol);

    std::string msg = "[NavierStokesAssemblerRB] computing rhs ...";
    printlog(YELLOW, msg, this->M_data.getVerbose());

    shp<BlockVector> retVec = convert<BlockVector>(StokesAssemblerRB::getRightHandSide(time,sol));

    bool approximatenonlinearterm = M_data("rb/online/approximatenonlinearterm", 1);

    if (!approximatenonlinearterm)
    {
        shp<VECTOREPETRA>  nonLinearTerm(new VECTOREPETRA(M_FEAssembler->getFEspace(0)->map()));
        shp<VECTOREPETRA>  velocityReconstructed;

        velocityReconstructed = M_bases->reconstructFEFunction(solBlck->block(0), 0, M_treeNode->M_ID);

        shp<ETFESPACE3> velocityFESpaceETA = M_FEAssembler->getVelocityETFEspace();
        double density = M_FEAssembler->getDensity();
        integrate(elements(velocityFESpaceETA->mesh()),
                   M_FEAssembler->getFEspace(0)->qr(),
                   velocityFESpaceETA,
                   value(density) *
                   dot(value(velocityFESpaceETA , *velocityReconstructed) *
                        grad(velocityFESpaceETA , *velocityReconstructed),
                   phi_i)
                 ) >> nonLinearTerm;

        nonLinearTerm->globalAssemble();

        shp<DistributedVector> nonLinearTermWrap(new DistributedVector());
        nonLinearTermWrap->setVector(nonLinearTerm);

        shp<BlockVector> nonLinearTermVec(new BlockVector(2));
        nonLinearTermVec->setBlock(0,nonLinearTermWrap);

        getBCManager()->apply0DirichletBCs(*nonLinearTermVec,
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
                shp<BlockVector> currVec(new BlockVector(0));
                currVec->deepCopy(M_nonLinearTermsDecomposition[i][j]);
                currVec->multiplyByScalar(solBlck->block(0)->operator()(i) *
                                          solBlck->block(0)->operator()(j));
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

shp<aMatrix>
NavierStokesAssemblerRB::
getJacobianRightHandSide(const double& time,
                         const shp<aVector>& sol)
{
    shp<aMatrix> jac = StokesAssemblerRB::getJacobianRightHandSide(time,sol);

    if (M_exactJacobian)
        addConvectiveTermJacobianRightHandSide(sol, nullptr, jac);

    return jac;
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

        auto velocityBasis = M_bases->getEnrichedBasis(0, M_treeNode->M_ID);

        int nterms = M_data("rb/online/numbernonlinearterms", 20);

        if (nterms == -1)
            nterms = velocityBasis.size();

        shp<MATRIXEPETRA> nonLinearMatrix(new MATRIXEPETRA(M_FEAssembler->getFEspace(0)->map()));
        shp<VECTOREPETRA> nonLinearTerm(new VECTOREPETRA(M_FEAssembler->getFEspace(0)->map()));

        M_nonLinearTermsDecomposition.resize(nterms);
        for (unsigned int i = 0; i < nterms; i++)
        {
            M_nonLinearTermsDecomposition[i].resize(nterms);
            nonLinearMatrix->zero();

            shp<ETFESPACE3> velocityFESpaceETA = M_FEAssembler->getVelocityETFEspace();
            double density = M_FEAssembler->getDensity();
            integrate(elements(velocityFESpaceETA->mesh()),
                       M_FEAssembler->getFEspace(0)->qr(),
                       velocityFESpaceETA,
                       velocityFESpaceETA,
                       value(density) *
                       dot(value(velocityFESpaceETA , *velocityBasis[i]) * grad(phi_j),
                       phi_i)
                     ) >> nonLinearMatrix;

            nonLinearMatrix->globalAssemble();

            for (unsigned int j = 0; j < nterms; j++)
            {
                *nonLinearTerm = (*nonLinearMatrix) * (*velocityBasis[j]);

                shp<DistributedVector> nonLinearTermWrap(new DistributedVector());
                nonLinearTermWrap->setVector(nonLinearTerm);

                shp<BlockVector> nonLinearTermVec(new BlockVector(2));
                nonLinearTermVec->setBlock(0,nonLinearTermWrap);

                getBCManager()->apply0DirichletBCs(*nonLinearTermVec,
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
