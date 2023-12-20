#include "NavierStokesAssemblerRB.hpp"

namespace RedMA
{

NavierStokesAssemblerRB::
NavierStokesAssemblerRB(const DataContainer& data,
                        shp<TreeNode> treeNode) :
  StokesAssemblerRB(data, treeNode)
{
    M_name = "NavierStokesAssemblerRB";
    M_exactJacobian = M_data("rb/online/exactjacobian",false);
    M_stabilizationName = data("assembler/stabilization/type", "none");

    // if we use a stabilization we use P1-P1 by default and exact jacobian
    if (std::strcmp(M_stabilizationName.c_str(), "none"))
    {
        M_FEAssembler->setVelocityOrder("P1");
        M_FEAssembler->setPressureOrder("P1");

        M_exactJacobian = true;
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

        std::string meshPath = M_treeNode->M_block->getMeshName();
        unsigned int dashpos = meshPath.find_last_of("/");
        unsigned int formatpos = meshPath.find_last_of(".");
        std::string nameMesh = meshPath.substr(dashpos + 1, meshPath.size()+dashpos+1-formatpos);

        auto velocityBasis = M_bases->getEnrichedBasis(0, ID());

        int nterms = M_data("rb/online/numbernonlinearterms", 10);

        if (nterms == -1)
            nterms = velocityBasis.size();

        shp<MATRIXEPETRA> nonLinearMatrix(new MATRIXEPETRA(M_FEAssembler->getFEspace(0)->map()));
        shp<MATRIXEPETRA> nonLinearJacobian(new MATRIXEPETRA(M_FEAssembler->getFEspace(0)->map()));
        shp<VECTOREPETRA> nonLinearTerm(new VECTOREPETRA(M_FEAssembler->getFEspace(0)->map()));

        auto velocityNorm = M_FEAssembler->getNorm(0, false);

        M_nonLinearTermsDecomposition.resize(nterms);
        M_nonLinearMatrixDecomposition.resize(nterms);

        std::string dir_mat = "NLterm/" + nameMesh + "/Matrix/Block" + std::to_string(ID());
        fs::create_directories(dir_mat);
        std::string dir_vec = "NLterm/" + nameMesh + "/Vector/Block" + std::to_string(ID());
        fs::create_directories(dir_vec);
        std::string dir_mat2 = "NLterm/" + nameMesh + "/Matrix_PG/Block" + std::to_string(ID());
        fs::create_directories(dir_mat2);
        std::string dir_vec2 = "NLterm/" + nameMesh + "/Vector_PG/Block" + std::to_string(ID());
        fs::create_directories(dir_vec2);

        double density = M_FEAssembler->getDensity();
        shp<ETFESPACE3 > velocityFESpaceETA = M_FEAssembler->getVelocityETFEspace();

        for (unsigned int i = 0; i < nterms; i++)
        {
            std::string msg =  "\nConsidering basis number " + std::to_string(i);
            printlog(WHITE, msg, this->M_data.getVerbose());

            M_nonLinearTermsDecomposition[i].resize(nterms);

            nonLinearJacobian->zero();
            integrate(elements(velocityFESpaceETA->mesh()),
                      M_FEAssembler->getFEspace(0)->qr(),
                      velocityFESpaceETA,
                      velocityFESpaceETA,
                      value(density) *
                      dot(value(velocityFESpaceETA, *velocityBasis[i]) * grad(phi_j) +
                      phi_j * grad(velocityFESpaceETA , *velocityBasis[i]),
                      phi_i)
            ) >> nonLinearJacobian;
            nonLinearJacobian->globalAssemble();

            shp<BlockMatrix> nonLinearJacobianVec(new BlockMatrix(2,2));
            nonLinearJacobianVec->setBlock(0,0,wrap(nonLinearJacobian));
            this->M_FEAssembler->applyDirichletBCsMatrix(nonLinearJacobianVec, 0.0);

            // TODO: I changed here!
            auto jac00 = M_bases->matrixProject(nonLinearJacobianVec->block(0,0), 0, 0, ID());
            // auto jac00 = M_bases->matrixProject(nonLinearJacobianVec->block(0,0), 0, 0, ID(), velocityNorm);
            auto jac00_2 = M_bases->rightProject(spcast<SparseMatrix>(nonLinearJacobianVec->block(0,0)), 0, ID());

            M_nonLinearMatrixDecomposition[i].reset(new BlockMatrix(2,2));
            M_nonLinearMatrixDecomposition[i]->setBlock(0,0,jac00);

            std::string filename_mat = dir_mat + "/Mat_" + std::to_string(i) + ".m";
            jac00->dump(filename_mat);
            std::string filename_mat2 = dir_mat2 + "/Mat_" + std::to_string(i) + ".m";
            jac00_2->dump(filename_mat2);

            nonLinearMatrix->zero();
            integrate(elements(velocityFESpaceETA->mesh()),
                      M_FEAssembler->getFEspace(0)->qr(),
                      velocityFESpaceETA,
                      velocityFESpaceETA,
                      value(density) *
                      dot(value(velocityFESpaceETA, *velocityBasis[i]) * grad(phi_j),
                      phi_i)
                      ) >> nonLinearMatrix;
            nonLinearMatrix->globalAssemble();

            for (unsigned int j = 0; j < nterms; j++)
            {
                *nonLinearTerm = (*nonLinearMatrix) * (*velocityBasis[j]);

                shp<BlockVector> nonLinearTermVec(new BlockVector(2));
                nonLinearTermVec->setBlock(0, wrap(nonLinearTerm));
                this->M_FEAssembler->apply0DirichletBCs(nonLinearTermVec);

                // TODO: I changed here!
                M_nonLinearTermsDecomposition[i][j].reset(new BlockVector(2));
                M_nonLinearTermsDecomposition[i][j]->setBlock(0, M_bases->leftProject(nonLinearTermVec->block(0), 0,  ID()));
                // M_nonLinearTermsDecomposition[i][j]->setBlock(0, M_bases->leftProject(nonLinearTermVec->block(0), 0, ID(), velocityNorm));

                std::string filename_vec = dir_vec + "/Vec_" + std::to_string(i) + "_" + std::to_string(j) + ".m";
                M_nonLinearTermsDecomposition[i][j]->block(0)->dump(filename_vec);
                std::string filename_vec2 = dir_vec2 + "/Vec_" + std::to_string(i) + "_" + std::to_string(j) + ".m";
                nonLinearTermVec->block(0)->dump(filename_vec2);
            }
        }

        std::string msg = "done, in ";
        msg += std::to_string(chrono.diff());
        msg += " seconds\n";
        printlog(YELLOW, msg, this->M_data.getVerbose());
    }

    if (!std::strcmp(M_stabilizationName.c_str(),"none"))
    {
        printlog(YELLOW, "[NavierStokesAssemblerRB] Proceeding without stabilization...\n",
                 this->M_data.getVerbose());
    }
    else
    {
        if (!std::strcmp(M_stabilizationName.c_str(), "supg"))
        {
            printlog(WHITE, "[NavierStokesAssemblerRB] Setting up SUPG stabilization...\n",
                     this->M_data.getVerbose());

            M_stabilization.reset(new SUPGStabilization(this->M_data,
                                                        this->M_FEAssembler->getFEspace(0),
                                                        this->M_FEAssembler->getFEspace(1),
                                                        this->M_FEAssembler->getVelocityETFEspace(),
                                                        this->M_FEAssembler->getPressureETFEspace(),
                                                        this->M_comm));
        }
        else if (!std::strcmp(M_stabilizationName.c_str(), "ip"))
        {
            printlog(YELLOW, "[NavierStokesAssemblerRB] Setting up IP stabilization...\n",
                     this->M_data.getVerbose());

            M_stabilization.reset(new IPStabilization(this->M_data,
                                                        this->M_FEAssembler->getFEspace(0),
                                                        this->M_FEAssembler->getFEspace(1),
                                                        this->M_FEAssembler->getVelocityETFEspace(),
                                                        this->M_FEAssembler->getPressureETFEspace(),
                                                        this->M_comm));
        }
        else
            throw new Exception("Stabilization " + M_stabilizationName + " is not implemented!");

        M_stabilization->setDensityAndViscosity(this->M_FEAssembler->getDensity(),
                                                this->M_FEAssembler->getViscosity());
        M_stabilization->setup();
    }
}

shp<aMatrix>
NavierStokesAssemblerRB::
getMass(const double& time,
        const shp<aVector>& sol)
{
    shp<BlockMatrix> retMat(new BlockMatrix(this->M_nComponents,this->M_nComponents));
    retMat->deepCopy(this->M_reducedMass);

// TODO: optimize with respect to the type of stabilization. With SUPG very few can be done, as it is
//  highly non-linear; a possible upgrade would consist in reconstructing the FEM solution from the
//  RB one just one for every Newton iteration and not every time a new matrix/vector has to be assembled.
//  With IP (if it works!), instead, a lot can be optimized; indeed the RB stabilization matrix can be
//  assembled once-for-all at offline phase and, during the online one, there's never need to reconstruct
//  FEM vectors from RB ones or to project FEM matrices/vectors onto the RB subspace.

    if (M_stabilization)
    {
        shp<BlockVector> solFEM = this->reconstructFESolution(convert<BlockVector>(sol));
        shp<BlockVector> rhsFEM = convert<BlockVector>(this->M_FEAssembler->getForcingTerm(time));

        shp<BlockMatrix> stabMatFEM = M_stabilization->getMass(solFEM, rhsFEM);
        this->M_FEAssembler->applyDirichletBCsMatrix(stabMatFEM, 0.0);

        auto stabMat00 = M_bases->matrixProject(stabMatFEM->block(0,0), 0, 0, ID());
        auto stabMat10 = M_bases->matrixProject(stabMatFEM->block(1,0), 1, 0, ID());

        retMat->block(0,0)->add(stabMat00);
        retMat->block(1,0)->add(stabMat10);
    }

    return retMat;
}

shp<aMatrix>
NavierStokesAssemblerRB::
getMassJacobian(const double& time,
                const shp<aVector>& sol)
{
// TODO: optimize with respect to the type of stabilization. With SUPG very few can be done, as it is
//  highly non-linear; a possible upgrade would consist in reconstructing the FEM solution from the
//  RB one just one for every Newton iteration and not every time a new matrix/vector has to be assembled.
//  With IP (if it works!), instead, a lot can be optimized; indeed the RB stabilization matrix can be
//  assembled once-for-all at offline phase and, during the online one, there's never need to reconstruct
//  FEM vectors from RB ones or to project FEM matrices/vectors onto the RB subspace.

    shp<BlockMatrix> retMat(new BlockMatrix(this->M_nComponents,
                                            this->M_nComponents));

    if (M_stabilization)
    {
        shp<BlockVector> solFEM = this->reconstructFESolution(convert<BlockVector>(sol));
        shp<BlockVector> rhsFEM = convert<BlockVector>(this->M_FEAssembler->getForcingTerm(time));

        shp<BlockMatrix> stabJacFEM = M_stabilization->getMassJacobian(solFEM, rhsFEM);
        this->M_FEAssembler->applyDirichletBCsMatrix(stabJacFEM, 0.0);

        auto stabJac00 = M_bases->matrixProject(stabJacFEM->block(0,0), 0, 0, ID());

        retMat->block(0,0)->add(stabJac00);
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

    using namespace LifeV;
    using namespace ExpressionAssembly;

    shp<VECTOREPETRA>  velocityHandler;
    bool solIsFEM = (3*M_FEAssembler->getFEspace(0)->dof().numTotalDof() ==
                     spcast<VECTOREPETRA>(convert<BlockVector>(sol)->block(0)->data())->size())
                     &&
                    (M_FEAssembler->getFEspace(1)->dof().numTotalDof() ==
                     spcast<VECTOREPETRA>(convert<BlockVector>(sol)->block(1)->data())->size());
    if (!solIsFEM)
        velocityHandler = M_bases->reconstructFEFunction(convert<BlockVector>(sol)->block(0), 0, ID());
    else
        velocityHandler.reset(new VECTOREPETRA(*spcast<VECTOREPETRA>(convert<BlockVector>(sol)->block(0)->data())));

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

    auto projectedMat = M_bases->matrixProject(wrap(convectiveMatrix), 0, 0, ID());

    projectedMat->multiplyByScalar(-1);
    convert<BlockMatrix>(mat)->block(0,0)->add(projectedMat);
}

shp<aVector>
NavierStokesAssemblerRB::
getRightHandSide(const double& time,
                 const shp<aVector>& sol)
{
// TODO: optimize with respect to the type of stabilization. With SUPG very few can be done, as it is
//  highly non-linear; a possible upgrade would consist in reconstructing the FEM solution from the
//  RB one just one for every Newton iteration and not every time a new matrix/vector has to be assembled.
//  With IP (if it works!), instead, a lot can be optimized; indeed the RB stabilization matrix can be
//  assembled once-for-all at offline phase and, during the online one, there's never need to reconstruct
//  FEM vectors from RB ones or to project FEM matrices/vectors onto the RB subspace.

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

        this->M_FEAssembler->apply0DirichletBCs(nonLinearTermVec);

        M_nonLinearTerm = M_bases->leftProject(nonLinearTermVec, ID());
    }
    else
    {
        unsigned int nterms = M_nonLinearTermsDecomposition.size();
        M_nonLinearTerm.reset(new BlockVector(2));
        M_nonLinearTerm->multiplyByScalar(0);
        for (unsigned int i = 0; i < nterms; i++)
        {
            for (unsigned int j = 0; j < nterms; j++)
            {
                shp<BlockVector> currVec(new BlockVector(0));
                currVec->deepCopy(M_nonLinearTermsDecomposition[i][j]);
                currVec->multiplyByScalar(solBlck->block(0)->operator()(i) *
                                          solBlck->block(0)->operator()(j));

                M_nonLinearTerm->block(0)->add(currVec->block(0));
            }
        }
    }

    M_nonLinearTerm->multiplyByScalar(-1);
    retVec->block(0)->add(M_nonLinearTerm->block(0));
    M_nonLinearTerm->multiplyByScalar(-1);

    if (M_stabilization)
    {
        shp<BlockVector> solFEM = this->reconstructFESolution(convert<BlockVector>(sol));
        shp<BlockVector> rhsFEM = convert<BlockVector>(this->M_FEAssembler->getForcingTerm(time));

        shp<BlockVector> residual = M_stabilization->getResidual(solFEM, rhsFEM);
        residual->multiplyByScalar(-1.0);

        residual = M_bases->leftProject(residual, ID());

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
// TODO: optimize with respect to the type of stabilization. With SUPG very few can be done, as it is
//  highly non-linear; a possible upgrade would consist in reconstructing the FEM solution from the
//  RB one just one for every Newton iteration and not every time a new matrix/vector has to be assembled.
//  With IP (if it works!), instead, a lot can be optimized; indeed the RB stabilization matrix can be
//  assembled once-for-all at offline phase and, during the online one, there's never need to reconstruct
//  FEM vectors from RB ones or to project FEM matrices/vectors onto the RB subspace.

    Chrono chrono;
    chrono.start();

    if (!std::strcmp((this->M_name).c_str(), "NavierStokesAssemblerRB"))
    {
        std::string msg = "[NavierStokesAssemblerRB] computing right-hand side jacobian ...";
        printlog(YELLOW, msg, this->M_data.getVerbose());
    }

    shp<aMatrix> jac = StokesAssemblerRB::getJacobianRightHandSide(time,sol);
    shp<BlockVector> solFEM;

    if (M_exactJacobian)
    {
        solFEM = this->reconstructFESolution(convert<BlockVector>(sol));
        addConvectiveTermJacobian(solFEM, jac);
    }

    if (M_stabilization)
    {
        if (!solFEM)
            solFEM = this->reconstructFESolution(convert<BlockVector>(sol));
        shp<BlockVector> rhsFEM = convert<BlockVector>(this->M_FEAssembler->getForcingTerm(time));

        shp<BlockMatrix> stabJacFEM = M_stabilization->getJacobian(solFEM, rhsFEM);
        stabJacFEM->multiplyByScalar(-1.0);
        this->M_FEAssembler->applyDirichletBCsMatrix(stabJacFEM, 0.0);

        auto stabJac00 = M_bases->matrixProject(stabJacFEM->block(0,0), 0, 0, ID());
        auto stabJac01 = M_bases->matrixProject(stabJacFEM->block(0,1), 0, 1, ID());
        auto stabJac10 = M_bases->matrixProject(stabJacFEM->block(1,0), 1, 0, ID());
        auto stabJac11 = M_bases->matrixProject(stabJacFEM->block(1,1), 1, 1, ID());

        jac->block(0,0)->add(stabJac00);
        jac->block(0,1)->add(stabJac01);
        jac->block(1,0)->add(stabJac10);
        jac->block(1,1)->add(stabJac11);
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
