#include "InterfaceAssembler.hpp"

namespace RedMA
{

template <>
std::vector<DenseVector>
InterfaceAssembler<DenseVector, DenseMatrix>::
buildCouplingVectors(SHP(BasisFunctionFunctor) bfs,
                     const GeometricFace& face,
                     SHP(aAssembler<DenseVector COMMA DenseMatrix>) assembler) const
{
    unsigned int nBasisFunctions = bfs->getNumBasisFunctions();
    std::vector<DenseVector> couplingVectors(3 * nBasisFunctions);

    return couplingVectors;
}

template <>
std::vector<DenseVector>
InterfaceAssembler<DenseVector, DenseMatrix>::
buildStabilizationVectorsVelocity(SHP(BasisFunctionFunctor) bfs,
                                  const GeometricFace& face,
                                  SHP(aAssembler<DenseVector COMMA DenseMatrix>) assembler) const
{
    unsigned int nBasisFunctions = bfs->getNumBasisFunctions();
    std::vector<DenseVector> couplingVectors(3 * nBasisFunctions);

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

    return couplingVectors;
}

template <>
void
InterfaceAssembler<DenseVector, DenseMatrix>::
buildCouplingMatrices(SHP(AssemblerType) assembler, const GeometricFace& face,
                      BlockMatrix<DenseMatrix>& matrixT, BlockMatrix<DenseMatrix>& matrix)
{
}

template <>
std::vector<DenseVector>
InterfaceAssembler<DenseVector, DenseMatrix>::
buildStabilizationVectorsLagrange() const
{
    std::vector<DenseVector> retVectors(0);

    return retVectors;
}

template <>
void
InterfaceAssembler<DenseVector, DenseMatrix>::
buildStabilizationMatrix(SHP(AssemblerType) assembler,
                         const GeometricFace& face,
                         BlockMatrix<DenseMatrix>& matrix)
{
}

template <>
void
InterfaceAssembler<DenseVector, DenseMatrix>::
buildCouplingMatrices()
{
}

template <>
BlockVector<DenseVector>
InterfaceAssembler<DenseVector, DenseMatrix>::
getZeroVector() const
{
    BlockVector<DenseVector> retVector;
    retVector.resize(1);

    return retVector;
}

}
