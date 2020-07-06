#include "InterfaceAssembler.hpp"

namespace RedMA
{

template <>
std::vector<DenseVector>
InterfaceAssembler<DenseVector, DenseMatrix>::
buildStabilizationVectorsVelocity(SHP(BasisFunctionFunctor) bfs,
                                  const GeometricFace& face,
                                  SHP(aAssembler<DenseVector COMMA DenseMatrix>) assembler) const
{
    unsigned int nBasisFunctions = bfs->getNumBasisFunctions();
    std::vector<DenseVector> couplingVectors(3 * nBasisFunctions);
    throw new Exception("Stabilization matrices not implemented in rb setting");

    return couplingVectors;
}

template <>
std::vector<DenseVector>
InterfaceAssembler<DenseVector, DenseMatrix>::
buildStabilizationVectorsPressure(SHP(BasisFunctionFunctor) bfs,
                                  const GeometricFace& face,
                                  SHP(aAssembler<DenseVector COMMA DenseMatrix>) assembler) const
{
    unsigned int nBasisFunctions = bfs->getNumBasisFunctions();
    std::vector<DenseVector> couplingVectors(3 * nBasisFunctions);

    throw new Exception("Stabilization matrices not implemented in rb setting");

    return couplingVectors;
}

template <>
std::vector<DenseVector>
InterfaceAssembler<DenseVector, DenseMatrix>::
buildStabilizationVectorsLagrange() const
{
    std::vector<DenseVector> retVectors(0);

    throw new Exception("Stabilization matrices not implemented in rb setting");

    return retVectors;
}

template <>
void
InterfaceAssembler<DenseVector, DenseMatrix>::
buildStabilizationMatrix(SHP(AssemblerType) assembler,
                         const GeometricFace& face,
                         BlockMatrix<DenseMatrix>& matrix)
{
    throw new Exception("Stabilization matrices not implemented in rb setting");
}

template <>
void
InterfaceAssembler<DenseVector, DenseMatrix>::
buildCouplingMatrices()
{
    unsigned int indexOutlet = M_interface.M_indexOutlet;

    auto asFather = M_interface.M_assemblerFather;

    if (asFather)
    {
        GeometricFace outlet = asFather->getTreeNode()->M_block->getOutlet(indexOutlet);

        BlockMatrix<MatrixEp> fatherBT;
        BlockMatrix<MatrixEp> fatherB;
        buildCouplingMatrices(asFather, outlet, fatherBT, fatherB);

        M_fatherBT = asFather->getRBBases()->leftProject(fatherBT, asFather->ID());
        M_fatherB = asFather->getRBBases()->rightProject(fatherB, asFather->ID());
    }

    auto asChild = M_interface.M_assemblerChild;
    if (asChild)
    {
        GeometricFace inlet = asChild->getTreeNode()->M_block->getInlet();
        // I invert the normal of the face such that it is the same as the outlet
        inlet.M_normal *= (-1.);

        BlockMatrix<MatrixEp> childBT;
        BlockMatrix<MatrixEp> childB;
        buildCouplingMatrices(asChild, inlet, childBT, childB);

        childB *= (-1.);
        childBT *= (-1.);

        M_childBEp = childB;

        M_childBT = asChild->getRBBases()->leftProject(childBT, asChild->ID());
        M_childB = asChild->getRBBases()->rightProject(childB, asChild->ID());
    }
}

template <>
BlockVector<DenseVector>
InterfaceAssembler<DenseVector, DenseMatrix>::
getZeroVector() const
{
    BlockVector<DenseVector> retVector;
    retVector.resize(1);

    unsigned int nlag;
    if (M_interface.M_assemblerFather)
        nlag = M_fatherBT.block(0,0).data()->N();
    else
        nlag = M_childBT.block(0,0).data()->N();

    retVector.block(0).data().reset(new DENSEVECTOR(nlag));
    retVector.block(0).data()->Scale(0.0);

    return retVector;
}

}
