#include "StokesAssembler.hpp"

namespace RedMA
{

template <>
void
StokesAssembler<DenseVector, DenseMatrix>::
apply0DirichletBCsMatrix(BlockMatrix<DenseMatrix>& matrix, double diagCoeff) const
{
}

template <>
void
StokesAssembler<DenseVector, DenseMatrix>::
apply0DirichletBCs(BlockVector<DenseVector>& vector) const
{
}

template <>
void
StokesAssembler<DenseVector, DenseMatrix>::
applyDirichletBCs(const double& time, BlockVector<DenseVector>& vector) const
{
}

template <>
BlockMatrix<DenseMatrix>
StokesAssembler<DenseVector,DenseMatrix>::
assembleStiffness(BlockMDEIMStructure* structure)
{
    BlockMatrix<DenseMatrix> stiffness;

    return stiffness;
}

template <>
BlockMatrix<DenseMatrix>
StokesAssembler<DenseVector,DenseMatrix>::
assembleMass(BlockMDEIMStructure* structure)
{
    BlockMatrix<DenseMatrix> mass;

    return mass;
}

template <>
BlockMatrix<DenseMatrix>
StokesAssembler<DenseVector,DenseMatrix>::
assembleDivergence(BlockMDEIMStructure* structure)
{
    BlockMatrix<DenseMatrix> divergence;

    return divergence;
}

template <>
void
StokesAssembler<DenseVector, DenseMatrix>::
exportSolution(const double& t, const BlockVector<DenseVector>& sol)
{
}

//
// template <>
// void
// StokesAssembler<VectorEp,MatrixEp>::
// addBackFlowStabilization(BlockVector<VectorEp> input,
//                          const BlockVector<VectorEp>& sol,
//                          const unsigned int& faceFlag)
// {
//     using namespace LifeV;
//     using namespace ExpressionAssembly;
//
//     SHP(VECTOREPETRA) vn(new VECTOREPETRA(*sol.block(0).data()));
//
//     *vn *= *M_flowRateVectors[faceFlag];
//
//     SHP(VECTOREPETRA) absvn(new VECTOREPETRA(*vn));
//     absvn->abs();
//
//     *vn -= *absvn;
//     *vn /= 2.0;
//
//     *vn *= *sol.block(0).data();
//
//     SHP(VECTOREPETRA) vnRepeated(new VECTOREPETRA(*vn, Repeated));
//     SHP(VECTOREPETRA) backflowStabRepeated(new VECTOREPETRA(vn->mapPtr(), Repeated));
//
//     QuadratureBoundary myBDQR(buildTetraBDQR(quadRuleTria7pt));
//
//     integrate(boundary(M_velocityFESpaceETA->mesh(), faceFlag),
//               myBDQR,
//               M_velocityFESpaceETA,
//               dot(value(M_velocityFESpaceETA, *vnRepeated), phi_i)
//           ) >> backflowStabRepeated;
//
//     backflowStabRepeated->globalAssemble();
//
//     *backflowStabRepeated *= 0.2 * M_density;
//
//     SHP(VECTOREPETRA) backflowStab(new VECTOREPETRA(*backflowStabRepeated, Unique));
//
//     *input.block(0).data() += *backflowStab;
// }
//
template <>
BlockVector<DenseVector>
StokesAssembler<DenseVector,DenseMatrix>::
getZeroVector() const
{
    BlockVector<DenseVector> retVec;

    return retVec;
}


//
// template <>
// SHP(MATRIXEPETRA)
// StokesAssembler<VectorEp, MatrixEp>::
// assembleFlowRateJacobian(const GeometricFace& face)
// {
//     using namespace LifeV;
//     using namespace ExpressionAssembly;
//
//     const double dropTolerance(2.0 * std::numeric_limits<double>::min());
//
//     SHP(MAPEPETRA) rangeMap = M_flowRateVectors[face.M_flag]->mapPtr();
//     EPETRACOMM comm = rangeMap->commPtr();
//
//
//     Epetra_Map epetraMap = M_flowRateVectors[face.M_flag]->epetraMap();
//     unsigned int numElements = epetraMap.NumMyElements();
//     unsigned int numGlobalElements = epetraMap.NumGlobalElements();
//
//     // this should be optimized
//     SHP(MATRIXEPETRA) flowRateJacobian;
//     flowRateJacobian.reset(new MATRIXEPETRA(M_velocityFESpace->map(), numGlobalElements, false));
//
//     // compute outer product of flowrate vector with itself
//     for (unsigned int j = 0; j < numGlobalElements; j++)
//     {
//         double myvaluecol = 0;
//
//         if (M_flowRateVectors[face.M_flag]->isGlobalIDPresent(j))
//             myvaluecol = M_flowRateVectors[face.M_flag]->operator[](j);
//
//         double valuecol = 0;
//         comm->SumAll(&myvaluecol, &valuecol, 1);
//
//         if (std::abs(valuecol) > dropTolerance)
//         {
//             for (unsigned int i = 0; i < numElements; i++)
//             {
//                 unsigned int gdof = epetraMap.GID(i);
//                 if (M_flowRateVectors[face.M_flag]->isGlobalIDPresent(gdof))
//                 {
//                     double valuerow = M_flowRateVectors[face.M_flag]->operator[](gdof);
//                     if (std::abs(valuerow * valuecol) > dropTolerance)
//                     {
//                         flowRateJacobian->addToCoefficient(gdof, j, valuerow * valuecol);
//                     }
//                 }
//             }
//         }
//
//     }
//
//     comm->Barrier();
//
//     flowRateJacobian->globalAssemble();
//
//     return flowRateJacobian;
// }
//
// template <>
// void
// StokesAssembler<VectorEp, MatrixEp>::
// assembleFlowRateJacobians()
// {
//     // assemble inflow flow rate vector
//     if (M_treeNode->isInletNode())
//     {
//         auto face = M_treeNode->M_block->getInlet();
//         M_flowRateJacobians[face.M_flag].resize(this->M_nComponents,
//                                                 this->M_nComponents);
//         M_flowRateJacobians[face.M_flag].block(0,0).data() = assembleFlowRateJacobian(face);
//         M_bcManager->apply0DirichletMatrix(M_flowRateJacobians[face.M_flag], getFESpaceBCs(),
//                                            getComponentBCs(), 0.0);
//     }
//
//     if (M_treeNode->isOutletNode())
//     {
//         auto faces = M_treeNode->M_block->getOutlets();
//
//         for (auto face : faces)
//         {
//             M_flowRateJacobians[face.M_flag].resize(this->M_nComponents,
//                                                     this->M_nComponents);
//             M_flowRateJacobians[face.M_flag].block(0,0).data() = assembleFlowRateJacobian(face);
//             M_bcManager->apply0DirichletMatrix(M_flowRateJacobians[face.M_flag], getFESpaceBCs(),
//                                                getComponentBCs(), 0.0);
//         }
//     }
// }
//
template <>
BlockVector<RBVECTOR>
StokesAssembler<RBVECTOR, RBMATRIX>::
getLifting(const double& time) const
{
    BlockVector<RBVECTOR> lifting;

    return lifting;
}

}
