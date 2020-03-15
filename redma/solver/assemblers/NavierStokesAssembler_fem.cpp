#include "NavierStokesAssembler.hpp"

namespace RedMA
{

template <>
void
NavierStokesAssembler<VectorEp,MatrixEp>::
addConvectiveMatrixRightHandSide(const BlockVector<VectorEp>& sol,
                                 BlockMatrix<MatrixEp>& mat)
{
    using namespace LifeV;
    using namespace ExpressionAssembly;

    SHP(MATRIXEPETRA)  convectiveMatrix(new MATRIXEPETRA(M_velocityFESpace->map()));
    SHP(VECTOREPETRA)  velocityRepeated(new VECTOREPETRA(*sol.block(0).data(),
                                                         Repeated));

    integrate(elements(M_velocityFESpaceETA->mesh()),
               M_velocityFESpace->qr(),
               M_velocityFESpaceETA,
               M_velocityFESpaceETA,
               value(this->M_density) *
               dot(value(M_velocityFESpaceETA , *velocityRepeated) * grad(phi_j),
               phi_i)
             ) >> convectiveMatrix;
    convectiveMatrix->globalAssemble();

    *mat.block(0,0).data() -= *convectiveMatrix;
}

template <>
void
NavierStokesAssembler<VectorEp,MatrixEp>::
addConvectiveTermJacobianRightHandSide(const BlockVector<VectorEp>& sol,
                                       const BlockVector<VectorEp>& lifting,
                                       BlockMatrix<MatrixEp>& mat)
{
    using namespace LifeV;
    using namespace ExpressionAssembly;

    SHP(MATRIXEPETRA)  convectiveMatrix(new MATRIXEPETRA(M_velocityFESpace->map()));
    SHP(VECTOREPETRA)  velocityRepeated(new VECTOREPETRA(*sol.block(0).data(),
                                                         Repeated));
    SHP(VECTOREPETRA)  liftingRepeated(new VECTOREPETRA(*lifting.block(0).data(),
                                                        Repeated));

    if (M_extrapolatedSolution.nRows() == 0)
    {
        integrate(elements(M_velocityFESpaceETA->mesh()),
                   M_velocityFESpace->qr(),
                   M_velocityFESpaceETA,
                   M_velocityFESpaceETA,
                   value(this->M_density) *
                   dot(
                   (
                   value(M_velocityFESpaceETA , *velocityRepeated) * grad(phi_j) +
                   phi_j * grad(M_velocityFESpaceETA , *velocityRepeated) +
                   value(M_velocityFESpaceETA , *liftingRepeated) * grad(phi_j)
                   ),
                   phi_i)
                 ) >> convectiveMatrix;
    }
    else
    {
        integrate(elements(M_velocityFESpaceETA->mesh()),
                   M_velocityFESpace->qr(),
                   M_velocityFESpaceETA,
                   M_velocityFESpaceETA,
                   value(this->M_density) *
                   dot(
                   (
                   value(M_velocityFESpaceETA , *velocityRepeated) * grad(phi_j) +
                   value(M_velocityFESpaceETA , *liftingRepeated) * grad(phi_j)
                   ),
                   phi_i)
                 ) >> convectiveMatrix;
    }
    convectiveMatrix->globalAssemble();

    *mat.block(0,0).data() -= *convectiveMatrix;
}

template <>
BlockMatrix<MatrixEp>
NavierStokesAssembler<VectorEp, MatrixEp>::
getMass(const double& time, const BlockVector<VectorEp>& sol)
{
    BlockMatrix<MatrixEp> retMat;
    retMat.hardCopy(this->M_mass);
    if (M_useStabilization)
    {
        retMat += M_stabilization->getMass(sol, this->getForcingTerm(time));
        this->M_bcManager->apply0DirichletMatrix(retMat, this->getFESpaceBCs(),
                                                 this->getComponentBCs(), 1.0);
    }

    return retMat;
}

template <>
BlockMatrix<MatrixEp>
NavierStokesAssembler<VectorEp, MatrixEp>::
getMassJacobian(const double& time, const BlockVector<VectorEp>& sol)
{
    BlockMatrix<MatrixEp> retMat(this->M_nComponents, this->M_nComponents);
    if (M_useStabilization)
    {
        retMat += M_stabilization->getMassJac(sol, this->getForcingTerm(time));
        // we do it here because matrices from stabilization have no bcs
        this->M_bcManager->apply0DirichletMatrix(retMat, this->getFESpaceBCs(),
                                                 this->getComponentBCs(), 0.0);
    }

    return retMat;
}

template <>
BlockVector<VectorEp>
NavierStokesAssembler<VectorEp, MatrixEp>::
getRightHandSide(const double& time, const BlockVector<VectorEp>& sol)
{
    BlockVector<VectorEp> retVec;
    BlockMatrix<MatrixEp> systemMatrix;

    systemMatrix.resize(this->M_nComponents, this->M_nComponents);
    systemMatrix += this->M_stiffness;
    systemMatrix += this->M_divergence;
    systemMatrix *= (-1.0);

    if (this->M_extrapolatedSolution.nRows() > 0)
        this->addConvectiveMatrixRightHandSide(this->M_extrapolatedSolution, systemMatrix);
    else
        this->addConvectiveMatrixRightHandSide(sol, systemMatrix);

    retVec.softCopy(systemMatrix * sol);

    this->addNeumannBCs(retVec, time, sol);

    if (M_useStabilization)
        retVec -= M_stabilization->getResidual(sol, this->getForcingTerm(time));

    return retVec;
}

template <>
BlockMatrix<MatrixEp>
NavierStokesAssembler<VectorEp, MatrixEp>::
getJacobianRightHandSide(const double& time,
                         const BlockVector<VectorEp>& sol)
{
    BlockMatrix<MatrixEp> retMat;
    retMat = StokesAssembler<VectorEp,MatrixEp>::getJacobianRightHandSide(time, sol);

    if (this->M_extrapolatedSolution.nRows() > 0)
        this->addConvectiveTermJacobianRightHandSide(this->M_extrapolatedSolution,
                                                     this->getZeroVector(), retMat);
    else
        this->addConvectiveTermJacobianRightHandSide(sol,
                                                     this->getZeroVector(), retMat);

    if (M_useStabilization)
        retMat -= M_stabilization->getJac(sol, this->getForcingTerm(time));

    this->M_bcManager->apply0DirichletMatrix(retMat, this->getFESpaceBCs(),
                                             this->getComponentBCs(), 0.0);
    return retMat;
}


}
