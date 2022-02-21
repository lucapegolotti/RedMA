#include "NavierStokesAssemblerFE.hpp"

namespace RedMA
{

NavierStokesAssemblerFE::
NavierStokesAssemblerFE(const DataContainer& data,
                        shp<TreeNode> treeNode) :
  StokesAssemblerFE(data,treeNode)
{
    M_stabilizationName = data("assembler/stabilization/type", "none");
    // if we use a stabilization we use P1-P1 by default
    if (std::strcmp(M_stabilizationName.c_str(), "none"))
    {
        setVelocityOrder("P1");
        setPressureOrder("P1");
    }
}

void
NavierStokesAssemblerFE::
setup()
{
    StokesAssemblerFE::setup();

    if (!std::strcmp(M_stabilizationName.c_str(),"none"))
    {
        printlog(YELLOW, "[NavierStokesAssemblerFE] Proceeding without stabilization...\n",
                 this->M_data.getVerbose());
    }
    else
    {
        if (!std::strcmp(M_stabilizationName.c_str(), "supg"))
        {
            printlog(YELLOW, "[NavierStokesAssemblerFE] Setting up SUPG stabilization...\n",
                     this->M_data.getVerbose());

            M_stabilization.reset(new SUPGStabilization(this->M_data,
                                                        this->M_velocityFESpace,
                                                        this->M_pressureFESpace,
                                                        this->M_velocityFESpaceETA,
                                                        this->M_pressureFESpaceETA,
                                                        this->M_comm));
        }
        else if (!std::strcmp(M_stabilizationName.c_str(), "ip"))
        {
            printlog(YELLOW, "[NavierStokesAssemblerFE] Setting up IP stabilization...\n",
                     this->M_data.getVerbose());

            M_stabilization.reset(new IPStabilization(this->M_data,
                                                        this->M_velocityFESpace,
                                                        this->M_pressureFESpace,
                                                        this->M_velocityFESpaceETA,
                                                        this->M_pressureFESpaceETA,
                                                        this->M_comm));
        }
        else
            throw new Exception("Stabilization " + M_stabilizationName + " is not implemented!");

        M_stabilization->setDensityAndViscosity(this->M_density, this->M_viscosity);
        M_stabilization->setup();
    }

}

void
NavierStokesAssemblerFE::
addConvectiveMatrix(shp<aVector> sol,
                    shp<aMatrix> mat)
{
    using namespace LifeV;
    using namespace ExpressionAssembly;

    shp<MATRIXEPETRA> convectiveMatrix(new MATRIXEPETRA(M_velocityFESpace->map()));
    shp<VECTOREPETRA> velocityHandler = spcast<VECTOREPETRA>(convert<BlockVector>(sol)->block(0)->data());
    shp<VECTOREPETRA> velocityRepeated(new VECTOREPETRA(*velocityHandler, Repeated));

    integrate(elements(M_velocityFESpaceETA->mesh()),
               M_velocityFESpace->qr(),
               M_velocityFESpaceETA,
               M_velocityFESpaceETA,
               value(this->M_density) *
               dot(value(M_velocityFESpaceETA , *velocityRepeated) * grad(phi_j),
               phi_i)
             ) >> convectiveMatrix;
    convectiveMatrix->globalAssemble();

    // *convectiveMatrix *= 0.0;

    *spcast<MATRIXEPETRA>(convert<BlockMatrix>(mat)->block(0,0)->data()) -= *convectiveMatrix;
}

void
NavierStokesAssemblerFE::
addConvectiveTermJacobian(shp<aVector> sol,
                          shp<aMatrix> mat)
{
    using namespace LifeV;
    using namespace ExpressionAssembly;

    shp<VECTOREPETRA> velocityHandler = spcast<VECTOREPETRA>(convert<BlockVector>(sol)->block(0)->data());

    shp<MATRIXEPETRA>  convectiveMatrix(new MATRIXEPETRA(M_velocityFESpace->map()));
    shp<VECTOREPETRA>  velocityRepeated(new VECTOREPETRA(*velocityHandler, Repeated));

    integrate(elements(M_velocityFESpaceETA->mesh()),
               M_velocityFESpace->qr(),
               M_velocityFESpaceETA,
               M_velocityFESpaceETA,
               value(this->M_density) *
               dot(
               (
               value(M_velocityFESpaceETA , *velocityRepeated) * grad(phi_j) +
               phi_j * grad(M_velocityFESpaceETA , *velocityRepeated)
               ),
               phi_i)
             ) >> convectiveMatrix;
    convectiveMatrix->globalAssemble();

    // *convectiveMatrix *= 0.0;

    *spcast<MATRIXEPETRA>(convert<BlockMatrix>(mat)->block(0,0)->data()) -= *convectiveMatrix;
}

shp<aMatrix>
NavierStokesAssemblerFE::
getMass(const double& time,
        const shp<aVector>& sol)
{
    shp<BlockMatrix> retMat(new BlockMatrix(this->M_nComponents,this->M_nComponents));
    retMat->deepCopy(this->M_mass);

    if (M_stabilization)
    {
        retMat->add(M_stabilization->getMass(convert<BlockVector>(sol),
                                             convert<BlockVector>(this->getForcingTerm(time))));

        this->M_bcManager->apply0DirichletMatrix(*retMat, this->getFESpaceBCs(),
                                                 this->getComponentBCs(), 1.0,
                                                 !(this->M_addNoSlipBC));
    }

    return retMat;
}

shp<aMatrix>
NavierStokesAssemblerFE::
getMassJacobian(const double& time,
                const shp<aVector>& sol)
{
    shp<BlockMatrix> retMat(new BlockMatrix(this->M_nComponents,
                                            this->M_nComponents));

    if (M_stabilization)
    {
        retMat->add(M_stabilization->getMassJacobian(spcast<BlockVector>(sol),
                                                     spcast<BlockVector>(this->getForcingTerm(time))));
        // we do it here because matrices from stabilization have no bcs
        this->M_bcManager->apply0DirichletMatrix(*retMat, this->getFESpaceBCs(),
                                                 this->getComponentBCs(), 0.0,
                                                 !(this->M_addNoSlipBC));
    }

    return retMat;
}

shp<aVector>
NavierStokesAssemblerFE::
getRightHandSide(const double& time,
                 const shp<aVector>& sol)
{
    shp<BlockMatrix> systemMatrix(new BlockMatrix(this->M_nComponents,
                                                    this->M_nComponents));

    systemMatrix->add(this->M_stiffness);
    systemMatrix->add(this->M_divergence);
    systemMatrix->multiplyByScalar(-1.0);

    this->addConvectiveMatrix(sol, systemMatrix);

    shp<aVector> retVec = systemMatrix->multiplyByVector(sol);

    if (M_stabilization)
    {
        shp<BlockVector> residual = M_stabilization->getResidual(spcast<BlockVector>(sol),
                                                                 spcast<BlockVector>(this->getForcingTerm(time)));
        residual->multiplyByScalar(-1.);
        retVec->add(residual);
    }

    StokesAssemblerFE::addNeumannBCs(time, sol, retVec);
    this->M_bcManager->apply0DirichletBCs(*spcast<BlockVector>(retVec), this->getFESpaceBCs(),
                                          this->getComponentBCs(), !(this->M_addNoSlipBC));

    return retVec;
}

shp<aMatrix>
NavierStokesAssemblerFE::
getJacobianRightHandSide(const double& time,
                         const shp<aVector>& sol,
                         const double& diagCoeff)
{
    shp<aMatrix> retMat = StokesAssemblerFE::getJacobianRightHandSide(time, sol);

    this->addConvectiveTermJacobian(sol, retMat);

    if (M_stabilization)
    {
        shp<BlockMatrix> stabJac = M_stabilization->getJacobian(spcast<BlockVector>(sol),
                                                                spcast<BlockVector>(this->getForcingTerm(time)));
        stabJac->multiplyByScalar(-1.);
        retMat->add(stabJac);
    }

    this->M_bcManager->apply0DirichletMatrix(*spcast<BlockMatrix>(retMat),
                                             this->getFESpaceBCs(),
                                             this->getComponentBCs(), diagCoeff,
                                             !(this->M_addNoSlipBC));

    return retMat;
}

}
