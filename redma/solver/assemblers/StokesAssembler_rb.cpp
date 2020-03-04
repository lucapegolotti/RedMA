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
assembleMass(BlockMDEIMStructure* structure)
{
    // load MDEIM
    // M_mdeimMass.reset(new BlockMDEIM(M_data, M_comm));


    BlockMatrix<DenseMatrix> mass;

    return mass;
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

template <>
BlockVector<DenseVector>
StokesAssembler<DenseVector,DenseMatrix>::
getZeroVector() const
{
    BlockVector<DenseVector> retVec;

    return retVec;
}

template <>
BlockVector<RBVECTOR>
StokesAssembler<RBVECTOR, RBMATRIX>::
getLifting(const double& time) const
{
    BlockVector<RBVECTOR> lifting;

    return lifting;
}

}
