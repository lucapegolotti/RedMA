#include "NavierStokesAssemblerRB.hpp"

namespace RedMA
{

NavierStokesAssemblerRB::
NavierStokesAssemblerRB(const DataContainer& data,
                        shp<TreeNode> treeNode,
                        std::string stabilizationName) :
  StokesAssemblerRB(data,treeNode),
  M_stabilizationName(stabilizationName)
{
    M_name = "NavierStokesAssemblerRB";
    M_exactJacobian = M_data("rb/online/exactjacobian",0);

    // if we use a stabilization we use P1-P1 by default
    if (std::strcmp(M_stabilizationName.c_str(),""))
    {
        M_FEAssembler->setVelocityOrder("P1");
        M_FEAssembler->setPressureOrder("P1");
    }
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
        printlog(YELLOW, "[NavierStokesAssemblerRB] precomputing non linear terms \t", M_data.getVerbose());
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

                shp<BlockVector> nonLinearTermVec(new BlockVector(2));
                nonLinearTermVec->setBlock(0, wrap(nonLinearTerm));

                getBCManager()->apply0DirichletBCs(*nonLinearTermVec,
                                                   getFESpaceBCs(),
                                                   getComponentBCs(),
                                                   !(this->M_FEAssembler->hasNoSlipBCs()));
                M_nonLinearTermsDecomposition[i][j] = M_bases->leftProject(nonLinearTermVec, ID());
            }
            M_nonLinearTerm.reset(new BlockVector(0));
        }

        if (!std::strcmp(M_stabilizationName.c_str(), "supg"))
        {
            printlog(WHITE, "[NavierStokesAssemblerRB] Setting up SUPG stabilization...\n",
                     this->M_data.getVerbose());

            M_stabilization.reset(new SUPGStabilization(this->M_data,
                                                        this->M_FEAssembler->getFEspace(0),
                                                        this->M_FEAssembler->getFEspace(1),
                                                        this->M_FEAssembler->getVelocityETFEspace(),
                                                        this->M_FEAssembler->getPressureETFEspace()));
            M_stabilization->setDensityAndViscosity(this->M_FEAssembler->getDensity(),
                                                    this->M_FEAssembler->getViscosity());
        }
        else if (!std::strcmp(M_stabilizationName.c_str(),""))
        {
            printlog(WHITE, "[NavierStokesAssemblerRB] Proceeding without stabilization...\n",
                     this->M_data.getVerbose());
        }
        else
            throw new Exception("Stabilization " + M_stabilizationName + " is not implemented!");

        std::string msg = "done, in ";
        msg += std::to_string(chrono.diff());
        msg += " seconds\n";
        printlog(YELLOW, msg, this->M_data.getVerbose());
    }
}

shp<aMatrix>
NavierStokesAssemblerRB::
getMass(const double& time,
        const shp<aVector>& sol)
{
    shp<BlockMatrix> retMat(new BlockMatrix(this->M_nComponents,this->M_nComponents));
    retMat->deepCopy(this->M_reducedMass);

// TODO: here the RB stabilization is implemented in a naive way; indeed every contribution
//  due to stabilization is computed and assembled in the FEM domain and then projected!
//  As basic step, it will be beneficial to store the FEM-reconstructed solution and rhs term
//  in the class, so that they do not have to be recomputed when the mass, mass jacobian, rhs and
//  rhs jacobian are assembled. However, it is necessary to store multiple terms (because of the
//  time-marching) and to update them at every Netwon iteration.

    if (M_stabilization)
    {
        shp<BlockVector> solFEM = this->reconstructFESolution(convert<BlockVector>(sol));
        shp<BlockVector> rhsFEM = this->reconstructFESolution(convert<BlockVector>(this->M_FEAssembler->getForcingTerm(time)));

        shp<BlockMatrix> stabMat = M_stabilization->getMass(solFEM, rhsFEM);

        this->M_bcManager->apply0DirichletMatrix(*stabMat, this->getFESpaceBCs(),
                                                 this->getComponentBCs(), 0.0,
                                                 !(this->M_FEAssembler->hasNoSlipBCs()));

        retMat->add(M_bases->matrixProject(stabMat, 0, 0, ID()));
    }

    return retMat;
}

shp<aMatrix>
NavierStokesAssemblerRB::
getMassJacobian(const double& time,
                const shp<aVector>& sol)
{
    shp<BlockMatrix> retMat(new BlockMatrix(this->M_nComponents,
                                            this->M_nComponents));

    if (M_stabilization)
    {
        shp<BlockVector> solFEM = this->reconstructFESolution(convert<BlockVector>(sol));
        shp<BlockVector> rhsFEM = this->reconstructFESolution(convert<BlockVector>(this->M_FEAssembler->getForcingTerm(time)));

        shp<BlockMatrix> stabMat = M_stabilization->getMassJacobian(solFEM, rhsFEM);

        this->M_bcManager->apply0DirichletMatrix(*stabMat, this->getFESpaceBCs(),
                                                 this->getComponentBCs(), 0.0,
                                                 !(this->M_FEAssembler->hasNoSlipBCs()));

        retMat->add(M_bases->matrixProject(stabMat, 0, 0, ID()));
    }

    return retMat;
}

void
NavierStokesAssemblerRB::
addConvectiveTermJacobian(shp<aVector> sol,
                          shp<aMatrix> mat)
{
    Chrono chrono;
    chrono.start();

    std::string msg = "[NavierStokesAssemblerRB] computing Jacobian rhs ...";
    printlog(YELLOW, msg, this->M_data.getVerbose());

    using namespace LifeV;
    using namespace ExpressionAssembly;

    shp<VECTOREPETRA>  velocityHandler;
    velocityHandler = M_bases->reconstructFEFunction(convert<BlockVector>(sol)->block(0), 0, ID());

    shp<MATRIXEPETRA>  convectiveMatrix(new MATRIXEPETRA(M_FEAssembler->getFEspace(0)->map()));
    shp<VECTOREPETRA>  velocityRepeated(new VECTOREPETRA(*velocityHandler, Repeated));

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

    unsigned int id = M_treeNode->M_ID;
    auto projectedMat = M_bases->matrixProject(wrap(convectiveMatrix), 0, 0, id);

    projectedMat->multiplyByScalar(-1);
    convert<BlockMatrix>(mat)->block(0,0)->add(projectedMat);

    msg = "done, in ";
    msg += std::to_string(chrono.diff());
    msg += " seconds\n";
    printlog(YELLOW, msg, this->M_data.getVerbose());
}

shp<aVector>
NavierStokesAssemblerRB::
getRightHandSide(const double& time,
                 const shp<aVector>& sol)
{
    using namespace LifeV;
    using namespace ExpressionAssembly;

    Chrono chrono;
    chrono.start();

    auto solBlck = convert<BlockVector>(sol);

    if (!std::strcmp((this->M_name).c_str(), "NavierStokesAssemblerRB"))
    {
        std::string msg = "[NavierStokesAssemblerRB] computing right-hand side term ...";
        printlog(YELLOW, msg, this->M_data.getVerbose());
    }

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

        shp<BlockVector> nonLinearTermVec(new BlockVector(2));
        nonLinearTermVec->setBlock(0, wrap(nonLinearTerm));

        getBCManager()->apply0DirichletBCs(*nonLinearTermVec,
                                           getFESpaceBCs(),
                                           getComponentBCs(),
                                           !(this->M_FEAssembler->hasNoSlipBCs()));

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
    retVec->block(0)->add(M_nonLinearTerm->block(0));
    M_nonLinearTerm->multiplyByScalar(-1);

    if (M_stabilization)
    {
        shp<BlockVector> solFEM = this->reconstructFESolution(convert<BlockVector>(sol));
        shp<BlockVector> rhsFEM = this->reconstructFESolution(convert<BlockVector>(this->M_FEAssembler->getForcingTerm(time)));

        shp<BlockVector> residual = M_stabilization->getResidual(solFEM, rhsFEM);
        residual->multiplyByScalar(-1.);
        retVec->add(residual);
    }

    if (!std::strcmp((this->M_name).c_str(), "NavierStokesAssemblerRB"))
    {
        std::string msg = "done, in ";
        msg += std::to_string(chrono.diff());
        msg += " seconds\n";
        printlog(YELLOW, msg, this->M_data.getVerbose());
    }

    return retVec;
}

shp<aMatrix>
NavierStokesAssemblerRB::
getJacobianRightHandSide(const double& time,
                         const shp<aVector>& sol)
{
    Chrono chrono;
    chrono.start();

    if (!std::strcmp((this->M_name).c_str(), "NavierStokesAssemblerRB"))
    {
        std::string msg = "[NavierStokesAssemblerRB] computing right-hand side jacobian ...";
        printlog(YELLOW, msg, this->M_data.getVerbose());
    }

    shp<aMatrix> jac = StokesAssemblerRB::getJacobianRightHandSide(time,sol);

    if (M_exactJacobian)
        addConvectiveTermJacobian(sol, jac);

    if (M_stabilization)
    {
        shp<BlockVector> solFEM = this->reconstructFESolution(convert<BlockVector>(sol));
        shp<BlockVector> rhsFEM = this->reconstructFESolution(convert<BlockVector>(this->M_FEAssembler->getForcingTerm(time)));

        shp<BlockMatrix> stabJac = M_stabilization->getJacobian(solFEM, rhsFEM);
        stabJac->multiplyByScalar(-1);
        jac->add(stabJac);
    }

    if (!std::strcmp((this->M_name).c_str(), "NavierStokesAssemblerRB"))
    {
        std::string msg = "done, in ";
        msg += std::to_string(chrono.diff());
        msg += " seconds\n";
        printlog(YELLOW, msg, this->M_data.getVerbose());
    }

    return jac;
}

}
