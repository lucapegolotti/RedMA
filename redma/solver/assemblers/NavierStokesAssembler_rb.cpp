#include "NavierStokesAssembler.hpp"

namespace RedMA
{

template <>
void
NavierStokesAssembler<DenseVector,DenseMatrix>::
addConvectiveMatrixRightHandSide(const BlockVector<DenseVector>& sol,
                                 BlockMatrix<DenseMatrix>& mat)
{
}

template <>
void
NavierStokesAssembler<DenseVector,DenseMatrix>::
addConvectiveTermJacobianRightHandSide(const BlockVector<DenseVector>& sol,
                                       const BlockVector<DenseVector>& lifting,
                                       BlockMatrix<DenseMatrix>& mat)
{
}

template <>
BlockMatrix<DenseMatrix>
NavierStokesAssembler<DenseVector, DenseMatrix>::
getMass(const double& time, const BlockVector<DenseVector>& sol)
{
    BlockMatrix<DenseMatrix> retMat;

    return retMat;
}

template <>
BlockMatrix<DenseMatrix>
NavierStokesAssembler<DenseVector, DenseMatrix>::
getMassJacobian(const double& time, const BlockVector<DenseVector>& sol)
{
    BlockMatrix<DenseMatrix> retMat(this->M_nComponents, this->M_nComponents);

    return retMat;
}

template <>
BlockVector<DenseVector>
NavierStokesAssembler<DenseVector, DenseMatrix>::
getRightHandSide(const double& time, const BlockVector<DenseVector>& sol)
{
    BlockVector<DenseVector> retVec;
    // BlockMatrix<InMatrixType> systemMatrix;
    //
    // systemMatrix.resize(this->M_nComponents, this->M_nComponents);
    // systemMatrix += this->M_stiffness;
    // systemMatrix += this->M_divergence;
    // systemMatrix *= (-1.0);
    //
    // if (this->extrapolatedSolution.nRows() > 0)
    //     this->addConvectiveMatrixRightHandSide(this->M_extrapolatedSolution, systemMatrix);
    // else
    //     this->addConvectiveMatrixRightHandSide(sol, systemMatrix);
    //
    // retVec.softCopy(systemMatrix * sol);
    //
    // this->addNeumannBCs(retVec, time, sol);
    //
    // if (M_useStabilization)
    //     retVec -= M_stabilization->getResidual(sol, this->getForcingTerm(time));

    return retVec;
}

template <>
BlockMatrix<DenseMatrix>
NavierStokesAssembler<DenseVector, DenseMatrix>::
getJacobianRightHandSide(const double& time,
                         const BlockVector<DenseVector>& sol)
{
    BlockMatrix<DenseMatrix> retMat;
    // retMat = StokesAssembler<InVectorType,InMatrixType>::getJacobianRightHandSide(time, sol);
    //
    // if (this->extrapolatedSolution.nRows() > 0)
    //     this->addConvectiveTermJacobianRightHandSide(this->M_extrapolatedSolution,
    //                                                  this->getZeroVector(), retMat);
    // else
    //     this->addConvectiveTermJacobianRightHandSide(sol,
    //                                                  this->getZeroVector(), retMat);
    //
    // if (M_useStabilization)
    //     retMat -= M_stabilization->getJac(sol, this->getForcingTerm(time));
    //
    // this->M_bcManager->apply0DirichletMatrix(retMat, this->getFESpaceBCs(),
    //                                          this->getComponentBCs(), 0.0);
    return retMat;
}


}
