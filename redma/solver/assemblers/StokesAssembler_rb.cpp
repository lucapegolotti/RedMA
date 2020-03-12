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
    M_mdeimMass->setFESpace(M_velocityFESpace, 0);
    M_mdeimMass->setFESpace(M_pressureFESpace, 1);

    BlockMatrix<MatrixEp> reducedMass = assembleReducedMass(&M_mdeimMass->getMDEIMStructure());
    BlockMatrix<DenseMatrix> mass = M_mdeimMass->assembleProjectedMatrix(reducedMass);

    return mass;
}

template <>
BlockMatrix<DenseMatrix>
StokesAssembler<DenseVector,DenseMatrix>::
assembleStiffness(BlockMDEIMStructure* structure)
{
    M_mdeimStiffness->setFESpace(M_velocityFESpace, 0);
    M_mdeimStiffness->setFESpace(M_pressureFESpace, 1);

    BlockMatrix<MatrixEp> reducedStiffness = assembleReducedStiffness(&M_mdeimStiffness->getMDEIMStructure());
    BlockMatrix<DenseMatrix> stiffness = M_mdeimStiffness->assembleProjectedMatrix(reducedStiffness);

    return stiffness;
}

template <>
BlockMatrix<DenseMatrix>
StokesAssembler<DenseVector,DenseMatrix>::
assembleDivergence(BlockMDEIMStructure* structure)
{
    M_mdeimDivergence->setFESpace(M_velocityFESpace, 0);
    M_mdeimDivergence->setFESpace(M_pressureFESpace, 1);

    BlockMatrix<MatrixEp> reducedDivergence = assembleReducedDivergence(&M_mdeimDivergence->getMDEIMStructure());
    BlockMatrix<DenseMatrix> divergence = M_mdeimDivergence->assembleProjectedMatrix(reducedDivergence);

    return divergence;
}

template <>
void
StokesAssembler<DenseVector, DenseMatrix>::
exportSolution(const double& t, const BlockVector<DenseVector>& sol)
{
    M_bases->reconstructFEFunction(M_velocityExporter, sol.block(0), 0);
    M_bases->reconstructFEFunction(M_pressureExporter, sol.block(1), 1);

    // BlockVector<VectorEp> solCopy(2);
    // solCopy.block(0).data() = M_velocityExporter;
    // computeFlowRates(solCopy, true);

    exportNorms(t);

    CoutRedirecter ct;
    ct.redirect();
    M_exporter->postProcess(t);
    printlog(CYAN, ct.restore());
}

template <>
BlockVector<DenseVector>
StokesAssembler<DenseVector,DenseMatrix>::
getZeroVector() const
{
    BlockVector<DenseVector> retVec;

    SHP(DENSEVECTOR) uComp(new DENSEVECTOR(M_bases->getSizeEnrichedBasis(0)));
    uComp->Scale(0.0);
    SHP(DENSEVECTOR) pComp(new DENSEVECTOR(M_bases->getSizeEnrichedBasis(1)));
    pComp->Scale(0.0);

    retVec.resize(M_nComponents);
    retVec.block(0).data() = uComp;
    retVec.block(1).data() = pComp;

    return retVec;
}

template <>
BlockVector<RBVECTOR>
StokesAssembler<RBVECTOR, RBMATRIX>::
getLifting(const double& time) const
{
    auto liftingFE = getFELifting(time);

    BlockVector<RBVECTOR> lifting;
    lifting = M_bases->leftProject(liftingFE);
    return lifting;
}

template <>
MatrixEp
StokesAssembler<DenseVector, DenseMatrix>::
getNorm(const unsigned int& fieldIndex)
{
    MatrixEp retMat;
    throw new Exception("Norm matrix not implemented for RB!");

    return retMat;
}

template <>
void
StokesAssembler<DenseVector, DenseMatrix>::
restrictRBMatrices()
{
    if (M_bases == nullptr)
        throw new Exception("RB bases have not been set yet");

    std::vector<unsigned int> selectorsU = M_bases->getSelectors(0);
    std::vector<unsigned int> selectorsP = M_bases->getSelectors(1);

    // this is the case when we do not choose to keep all the vectors in the basis
    if (selectorsU.size() > 0)
    {
        unsigned int Nu = selectorsU.size();
        unsigned int Np = selectorsP.size();

        // restrict mass
        SHP(DENSEMATRIX) restrictedMass(new DENSEMATRIX(Nu,Nu));

        unsigned int inew = 0;
        for (auto i : selectorsU)
        {
            unsigned int jnew = 0;
            for (auto j : selectorsU)
            {
                (*restrictedMass)(inew,jnew) = (*M_mass.block(0,0).data())(i,j);
                jnew++;
            }
            inew++;
        }

        M_mass.block(0,0).data() = restrictedMass;

        // restrict stiffness
        SHP(DENSEMATRIX) restrictedStiffness(new DENSEMATRIX(Nu,Nu));

        inew = 0;
        for (auto i : selectorsU)
        {
            unsigned int jnew = 0;
            for (auto j : selectorsU)
            {
                (*restrictedStiffness)(inew,jnew) = (*M_stiffness.block(0,0).data())(i,j);
                jnew++;
            }
            inew++;
        }

        M_stiffness.block(0,0).data() = restrictedStiffness;

        // restrict divergence
        SHP(DENSEMATRIX) restrictedBT(new DENSEMATRIX(Nu,Np));

        inew = 0;
        for (auto i : selectorsU)
        {
            unsigned int jnew = 0;
            for (auto j : selectorsP)
            {
                (*restrictedBT)(inew,jnew) = (*M_divergence.block(0,1).data())(i,j);
                jnew++;
            }
            inew++;
        }

        M_divergence.block(0,1).data() = restrictedBT;

        SHP(DENSEMATRIX) restrictedB(new DENSEMATRIX(Np,Nu));

        inew = 0;
        for (auto i : selectorsP)
        {
            unsigned int jnew = 0;
            for (auto j : selectorsU)
            {
                (*restrictedB)(inew,jnew) = (*M_divergence.block(1,0).data())(i,j);
                jnew++;
            }
            inew++;
        }

        M_divergence.block(1,0).data() = restrictedB;
    }
}

}
