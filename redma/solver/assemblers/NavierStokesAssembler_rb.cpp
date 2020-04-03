#include "NavierStokesAssembler.hpp"

namespace RedMA
{

template <>
void
NavierStokesAssembler<DenseVector,DenseMatrix>::
addConvectiveMatrixRightHandSide(const BlockVector<DenseVector>& sol,
                                 BlockMatrix<DenseMatrix>& mat)
{
    LifeV::LifeChrono chrono;
    chrono.start();

    std::string msg = "[NavierStokesAssembler] assembling convective matrix ...";
    printlog(YELLOW, msg, this->M_data.getVerbose());

    using namespace LifeV;
    using namespace ExpressionAssembly;

    SHP(MATRIXEPETRA)  convectiveMatrix(new MATRIXEPETRA(M_velocityFESpace->map()));
    SHP(VECTOREPETRA)  velocityRepeated;

    velocityRepeated = M_bases->reconstructFEFunction(sol.block(0), 0);

    integrate(elements(M_velocityFESpaceETA->mesh()),
               M_velocityFESpace->qr(),
               M_velocityFESpaceETA,
               M_velocityFESpaceETA,
               value(this->M_density) *
               dot(value(M_velocityFESpaceETA , *velocityRepeated) * grad(phi_j),
               phi_i)
             ) >> convectiveMatrix;
    convectiveMatrix->globalAssemble();

    BlockMatrix<MatrixEp> convectiveMatrixWrap(2,2);
    convectiveMatrixWrap.block(0,0).data() = convectiveMatrix;

    this->M_bcManager->apply0DirichletMatrix(convectiveMatrixWrap, this->getFESpaceBCs(),
                                             this->getComponentBCs(), 0.0);

    DenseMatrix convectiveMatrixProjected = M_bases->matrixProject(convectiveMatrixWrap.block(0,0),
                                                                   0, 0);

    mat.block(0,0) -= convectiveMatrixProjected;

    msg = "done, in ";
    msg += std::to_string(chrono.diff());
    msg += " seconds\n";
    printlog(YELLOW, msg, this->M_data.getVerbose());
}

template <>
void
NavierStokesAssembler<DenseVector,DenseMatrix>::
addConvectiveTermJacobianRightHandSide(const BlockVector<DenseVector>& sol,
                                       const BlockVector<DenseVector>& lifting,
                                       BlockMatrix<DenseMatrix>& mat)
{
    LifeV::LifeChrono chrono;
    chrono.start();

    std::string msg = "[NavierStokesAssembler] assembling convective jacobian matrix ...";
    printlog(YELLOW, msg, this->M_data.getVerbose());

    using namespace LifeV;
    using namespace ExpressionAssembly;

    SHP(MATRIXEPETRA)  convectiveMatrix(new MATRIXEPETRA(M_velocityFESpace->map()));
    SHP(VECTOREPETRA)  velocityRepeated;

    velocityRepeated = M_bases->reconstructFEFunction(sol.block(0), 0);

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

    BlockMatrix<MatrixEp> convectiveMatrixWrap(2,2);
    convectiveMatrixWrap.block(0,0).data() = convectiveMatrix;

    this->M_bcManager->apply0DirichletMatrix(convectiveMatrixWrap, this->getFESpaceBCs(),
                                             this->getComponentBCs(), 0.0);

    DenseMatrix convectiveMatrixProjected = M_bases->matrixProject(convectiveMatrixWrap.block(0,0), 0, 0);

    mat.block(0,0) -= convectiveMatrixProjected;

    msg = "done, in ";
    msg += std::to_string(chrono.diff());
    msg += " seconds\n";
    printlog(YELLOW, msg, this->M_data.getVerbose());
}

template <>
BlockMatrix<DenseMatrix>
NavierStokesAssembler<DenseVector, DenseMatrix>::
getMass(const double& time, const BlockVector<DenseVector>& sol)
{
    BlockMatrix<DenseMatrix> retMat;
    retMat.hardCopy(this->M_mass);

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

// template <>
// BlockVector<DenseVector>
// NavierStokesAssembler<DenseVector, DenseMatrix>::
// getRightHandSide(const double& time, const BlockVector<DenseVector>& sol)
// {
//     BlockVector<DenseVector> retVec;
//     BlockMatrix<DenseMatrix> systemMatrix;
//
//     systemMatrix.resize(this->M_nComponents, this->M_nComponents);
//     systemMatrix += this->M_stiffness;
//     systemMatrix += this->M_divergence;
//     systemMatrix *= (-1.0);
//
//     this->addConvectiveMatrixRightHandSide(sol, systemMatrix);
//
//     retVec.softCopy(systemMatrix * sol);
//
//     // this->addNeumannBCs(retVec, time, sol);
//
//     return retVec;
// }

template <>
BlockVector<DenseVector>
NavierStokesAssembler<DenseVector, DenseMatrix>::
getRightHandSide(const double& time, const BlockVector<DenseVector>& sol)
{
    LifeV::LifeChrono chrono;
    chrono.start();

    std::string msg = "[NavierStokesAssembler] computing rhs ...";
    printlog(YELLOW, msg, this->M_data.getVerbose());

    using namespace LifeV;
    using namespace ExpressionAssembly;

    BlockVector<DenseVector> retVec;
    BlockMatrix<DenseMatrix> systemMatrix;

    systemMatrix.resize(this->M_nComponents, this->M_nComponents);
    systemMatrix += this->M_stiffness;
    systemMatrix += this->M_divergence;
    systemMatrix *= (-1.0);

    retVec.softCopy(systemMatrix * sol);

    std::cout << systemMatrix.block(0,0).data()->M() << " " << systemMatrix.block(0,0).data()->N() << std::endl << std::flush;

    SHP(VECTOREPETRA)  nonLinearTerm(new VECTOREPETRA(M_velocityFESpace->map()));
    SHP(VECTOREPETRA)  velocityReconstructed;

    velocityReconstructed = M_bases->reconstructFEFunction(sol.block(0), 0);

    if (M_extrapolatedSolution.nRows() == 0)
    {
        integrate(elements(M_velocityFESpaceETA->mesh()),
                   M_velocityFESpace->qr(),
                   M_velocityFESpaceETA,
                   value(this->M_density) *
                   dot(value(M_velocityFESpaceETA , *velocityReconstructed) *
                       grad(M_velocityFESpaceETA , *velocityReconstructed),
                   phi_i)
                 ) >> nonLinearTerm;
    }
    else
    {
        SHP(VECTOREPETRA)  extrapolatedSolution;
        extrapolatedSolution = M_bases->reconstructFEFunction(M_extrapolatedSolution.block(0), 0);

        integrate(elements(M_velocityFESpaceETA->mesh()),
                   M_velocityFESpace->qr(),
                   M_velocityFESpaceETA,
                   value(this->M_density) *
                   dot(value(M_velocityFESpaceETA , *velocityReconstructed) *
                       grad(M_velocityFESpaceETA , *extrapolatedSolution),
                   phi_i)
                 ) >> nonLinearTerm;
    }

    nonLinearTerm->globalAssemble();

    BlockVector<VectorEp> nonLinearTermWrap(2);
    nonLinearTermWrap.block(0).data() = nonLinearTerm;

    if (M_useStabilization)
    {
        if (this->M_extrapolatedSolution.norm2() > 1e-15)
            throw new Exception("Stabilization is not supported with extrapolation");
        else
        {
            SHP(VECTOREPETRA)  pressureReconstructed;
            pressureReconstructed = M_bases->reconstructFEFunction(sol.block(1), 1);

            BlockVector<VectorEp> zeroVec(2);
            zeroVec.block(0).data().reset(new VECTOREPETRA(M_velocityFESpace->map()));
            zeroVec.block(1).data().reset(new VECTOREPETRA(M_pressureFESpace->map()));
            zeroVec *= 0.0;

            BlockVector<VectorEp> solWrap(2);
            solWrap.block(0).data() = velocityReconstructed;
            solWrap.block(1).data() = pressureReconstructed;

            nonLinearTermWrap += M_stabilization->getResidual(solWrap, zeroVec);
        }
    }

    M_bcManager->apply0DirichletBCs(nonLinearTermWrap,
                                    getFESpaceBCs(),
                                    getComponentBCs());

    BlockVector<DenseVector> nonLinearTermProjected = M_bases->leftProject(nonLinearTermWrap);


    retVec -= nonLinearTermProjected;
    // this->addNeumannBCs(retVec, time, sol);

    msg = "done, in ";
    msg += std::to_string(chrono.diff());
    msg += " seconds\n";
    printlog(YELLOW, msg, this->M_data.getVerbose());

    return retVec;
}

template <>
BlockMatrix<DenseMatrix>
NavierStokesAssembler<DenseVector, DenseMatrix>::
getJacobianRightHandSide(const double& time,
                         const BlockVector<DenseVector>& sol)
{
    BlockMatrix<DenseMatrix> retMat;
    retMat = StokesAssembler<DenseVector,DenseMatrix>::getJacobianRightHandSide(time, sol);

    // this->addConvectiveTermJacobianRightHandSide(sol, this->getZeroVector(), retMat);

    return retMat;
}


}
