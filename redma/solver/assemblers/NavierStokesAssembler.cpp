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
    M_bcManager->apply0DirichletMatrix(mat, getFESpaceBCs(),
                                       getComponentBCs(), 0.0);
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
    convectiveMatrix->globalAssemble();

    *mat.block(0,0).data() -= *convectiveMatrix;
    M_bcManager->apply0DirichletMatrix(mat, getFESpaceBCs(),
                                       getComponentBCs(), 0.0);
}

}
