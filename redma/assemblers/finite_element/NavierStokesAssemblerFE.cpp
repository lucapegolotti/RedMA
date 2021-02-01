#include "NavierStokesAssemblerFE.hpp"

namespace RedMA
{


NavierStokesAssemblerFE::
NavierStokesAssemblerFE(const DataContainer& data, shp<TreeNode> treeNode,
                        std::string stabilizationName) :
  StokesAssemblerFE(data,treeNode),
  NavierStokesModel(data,treeNode),
  M_stabilizationName(stabilizationName)
{
    // if we use a stabilization we use P1P1 by default
    if (std::strcmp(M_stabilizationName.c_str(),""))
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

    if (!std::strcmp(M_stabilizationName.c_str(),"supg"))
    {
        M_stabilization.reset(new SUPGStabilization(this->M_data,
                                                    this->M_velocityFESpace,
                                                    this->M_pressureFESpace,
                                                    this->M_velocityFESpaceETA,
                                                    this->M_pressureFESpaceETA));
        M_stabilization->setDensityAndViscosity(this->M_density, this->M_viscosity);
    }
    // else if (!std::strcmp(stabilizationType.c_str(),"vms"))
    // {
    //     M_stabilization.reset(new VMSStabilization(this->M_data,
    //                                                this->M_velocityFESpace,
    //                                                this->M_pressureFESpace,
    //                                                this->M_velocityFESpaceETA,
    //                                                this->M_pressureFESpaceETA));
    //     M_stabilization->setDensityAndViscosity(this->M_density, this->M_viscosity);
    // }
    // else if (!std::strcmp(stabilizationType.c_str(),"hf"))
    // {
    //     M_stabilization.reset(new HFStabilization(this->M_data,
    //                                               this->M_velocityFESpace,
    //                                               this->M_pressureFESpace,
    //                                               this->M_velocityFESpaceETA,
    //                                               this->M_pressureFESpaceETA));
    //     M_stabilization->setDensityAndViscosity(this->M_density, this->M_viscosity);
    // }
    else if (!std::strcmp(M_stabilizationName.c_str(),""))
    {

    }
    else
        throw new Exception("Type of stabilization not implemented!");
}

void
NavierStokesAssemblerFE::
addConvectiveMatrixRightHandSide(shp<aVector> sol, shp<aMatrix> mat)
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

    *spcast<MATRIXEPETRA>(convert<BlockMatrix>(mat)->block(0,0)->data()) -= *convectiveMatrix;
}

void
NavierStokesAssemblerFE::
addConvectiveTermJacobianRightHandSide(shp<aVector> sol, shp<aVector> lifting,
                                       shp<aMatrix> mat)
{
    using namespace LifeV;
    using namespace ExpressionAssembly;

    shp<VECTOREPETRA> velocityHandler = spcast<VECTOREPETRA>(convert<BlockVector>(sol)->block(0)->data());
    // shp<VECTOREPETRA> liftingHandler = spcast<VECTOREPETRA>(lifting->block(0)->data());

    shp<MATRIXEPETRA>  convectiveMatrix(new MATRIXEPETRA(M_velocityFESpace->map()));
    shp<VECTOREPETRA>  velocityRepeated(new VECTOREPETRA(*velocityHandler, Repeated));
    // shp<VECTOREPETRA>  liftingRepeated(new VECTOREPETRA(*liftingHandler, Repeated));

    // if the extrapolation is null (e.g. first step), the matrix is singular.
    // Hence we solve the non linear problem for the first step
    // if (M_extrapolatedSolution.norm2() > 1e-15)
    // {
    //     integrate(elements(M_velocityFESpaceETA->mesh()),
    //                M_velocityFESpace->qr(),
    //                M_velocityFESpaceETA,
    //                M_velocityFESpaceETA,
    //                value(this->M_density) *
    //                dot(
    //                (
    //                value(M_velocityFESpaceETA , *velocityRepeated) * grad(phi_j) +
    //                value(M_velocityFESpaceETA , *liftingRepeated) * grad(phi_j)
    //                ),
    //                phi_i)
    //              ) >> convectiveMatrix;
    // }
    // else
    // {
        integrate(elements(M_velocityFESpaceETA->mesh()),
                   M_velocityFESpace->qr(),
                   M_velocityFESpaceETA,
                   M_velocityFESpaceETA,
                   value(this->M_density) *
                   dot(
                   (
                   value(M_velocityFESpaceETA , *velocityRepeated) * grad(phi_j) +
                   phi_j * grad(M_velocityFESpaceETA , *velocityRepeated)
                   // +
                   // value(M_velocityFESpaceETA , *liftingRepeated) * grad(phi_j)
                   ),
                   phi_i)
                 ) >> convectiveMatrix;
    // }
    convectiveMatrix->globalAssemble();
    *spcast<MATRIXEPETRA>(convert<BlockMatrix>(mat)->block(0,0)->data()) -= *convectiveMatrix;
}

shp<aMatrix>
NavierStokesAssemblerFE::
getMass(const double& time, const shp<aVector>& sol)
{
    shp<BlockMatrix> retMat(new BlockMatrix(0,0));
    retMat->deepCopy(this->M_mass);
    if (M_stabilization)
    {
        // if (M_extrapolatedSolution.norm2() < 1e-15)
        retMat->add(M_stabilization->getMass(convert<BlockVector>(sol),
                                             convert<BlockVector>(this->getForcingTerm(time))));
        // else
        //retMat += M_stabilization->getMass(M_extrapolatedSolution, this->getForcingTerm(time));

        this->M_bcManager->apply0DirichletMatrix(*retMat, this->getFESpaceBCs(),
                                                 this->getComponentBCs(), 1.0,
                                                 !(this->M_addNoSlipBC));
    }

    return retMat;
}

shp<aMatrix>
NavierStokesAssemblerFE::
getMassJacobian(const double& time, const shp<aVector>& sol)
{
    shp<BlockMatrix> retMat(new BlockMatrix(this->M_nComponents,this->M_nComponents));
    if (M_stabilization)
    {
        retMat->add(M_stabilization->getMassJac(spcast<BlockVector>(sol),
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
getRightHandSide(const double& time, const shp<aVector>& sol)
{
    shp<BlockMatrix> systemMatrix(new BlockMatrix(this->M_nComponents,
                                                    this->M_nComponents));

    systemMatrix->add(M_stiffness);
    systemMatrix->add(M_divergence);
    systemMatrix->multiplyByScalar(-1.0);

    // if (this->M_extrapolatedSolution.nRows() > 0 && this->M_extrapolatedSolution.norm2() > 1e-15)
    //     this->addConvectiveMatrixRightHandSide(this->M_extrapolatedSolution, systemMatrix);
    // else
    this->addConvectiveMatrixRightHandSide(sol, systemMatrix);

    shp<aVector> retVec = systemMatrix->multiplyByVector(sol);

    addNeumannBCs(time, sol, retVec);
    this->M_bcManager->apply0DirichletBCs(*spcast<BlockVector>(retVec), this->getFESpaceBCs(),
                                          this->getComponentBCs(), !(this->M_addNoSlipBC));

    if (M_stabilization)
    {
        // if (this->M_extrapolatedSolution.norm2() > 1e-15)
        // {
        //     throw new Exception("Stabilization is not supported with extrapolation");
        //
        //     retVec -= M_stabilization->getResidual(M_extrapolatedSolution,
        //                                            this->getForcingTerm(time));
        // }
        // else
        shp<BlockVector> residual = M_stabilization->getResidual(spcast<BlockVector>(sol),
                                                                 spcast<BlockVector>(this->getForcingTerm(time)));
        residual->multiplyByScalar(-1);
        retVec->add(residual);
    }

    return retVec;
}

shp<aMatrix>
NavierStokesAssemblerFE::
getJacobianRightHandSide(const double& time,
                         const shp<aVector>& sol)
{
    shp<aMatrix> retMat = StokesAssemblerFE::getJacobianRightHandSide(time, sol);

    // if (this->M_extrapolatedSolution.nRows() > 0 && this->M_extrapolatedSolution.norm2() > 1e-15)
    //     this->addConvectiveTermJacobianRightHandSide(this->M_extrapolatedSolution,
    //                                                  this->getZeroVector(), retMat);
    // else
    this->addConvectiveTermJacobianRightHandSide(sol, this->getZeroVector(), retMat);

    if (M_stabilization)
    {
        shp<BlockMatrix> stabJac = M_stabilization->getJac(spcast<BlockVector>(sol),
                                                           spcast<BlockVector>(this->getForcingTerm(time)));
        stabJac->multiplyByScalar(-1);
        retMat->add(stabJac);
    }

    this->M_bcManager->apply0DirichletMatrix(*spcast<BlockMatrix>(retMat),
                                             this->getFESpaceBCs(),
                                             this->getComponentBCs(), 0.0,
                                             !(this->M_addNoSlipBC));

    return retMat;
}

}
