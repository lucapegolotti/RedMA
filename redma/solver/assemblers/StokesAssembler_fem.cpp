#include "StokesAssembler.hpp"

namespace RedMA
{

template <>
void
StokesAssembler<VectorEp, MatrixEp>::
apply0DirichletBCsMatrix(BlockMatrix<MatrixEp>& matrix, double diagCoeff) const
{
    this->M_bcManager->apply0DirichletMatrix(matrix, getFESpaceBCs(),
                                             getComponentBCs(), diagCoeff);
}

template <>
void
StokesAssembler<VectorEp, MatrixEp>::
apply0DirichletBCs(BlockVector<VectorEp>& vector) const
{
    this->M_bcManager->apply0DirichletBCs(vector, getFESpaceBCs(), getComponentBCs());
}

template <>
void
StokesAssembler<VectorEp, MatrixEp>::
applyDirichletBCs(const double& time, BlockVector<VectorEp>& vector) const
{
    this->M_bcManager->applyDirichletBCs(time, vector, getFESpaceBCs(),
                                         getComponentBCs());
}

template <>
BlockMatrix<MatrixEp>
StokesAssembler<VectorEp,MatrixEp>::
assembleStiffness(BlockMDEIMStructure* structure)
{
    return assembleReducedStiffness(structure);
}

template <>
BlockMatrix<MatrixEp>
StokesAssembler<VectorEp,MatrixEp>::
assembleMass(BlockMDEIMStructure* structure)
{
    return assembleReducedMass(structure);
}

template <>
BlockMatrix<MatrixEp>
StokesAssembler<VectorEp,MatrixEp>::
assembleDivergence(BlockMDEIMStructure* structure)
{
    return assembleReducedDivergence(structure);
}

template <>
void
StokesAssembler<VectorEp,MatrixEp>::
addBackFlowStabilization(BlockVector<VectorEp> input,
                         const BlockVector<VectorEp>& sol,
                         const unsigned int& faceFlag)
{
    using namespace LifeV;
    using namespace ExpressionAssembly;

    SHP(VECTOREPETRA) vn(new VECTOREPETRA(*sol.block(0).data()));

    *vn *= *M_flowRateVectors[faceFlag];

    SHP(VECTOREPETRA) absvn(new VECTOREPETRA(*vn));
    absvn->abs();

    *vn -= *absvn;
    *vn /= 2.0;

    *vn *= *sol.block(0).data();

    SHP(VECTOREPETRA) vnRepeated(new VECTOREPETRA(*vn, Repeated));
    SHP(VECTOREPETRA) backflowStabRepeated(new VECTOREPETRA(vn->mapPtr(), Repeated));

    QuadratureBoundary myBDQR(buildTetraBDQR(quadRuleTria7pt));

    integrate(boundary(M_velocityFESpaceETA->mesh(), faceFlag),
              myBDQR,
              M_velocityFESpaceETA,
              dot(value(M_velocityFESpaceETA, *vnRepeated), phi_i)
          ) >> backflowStabRepeated;

    backflowStabRepeated->globalAssemble();

    *backflowStabRepeated *= 0.2 * M_density;

    SHP(VECTOREPETRA) backflowStab(new VECTOREPETRA(*backflowStabRepeated, Unique));

    *input.block(0).data() += *backflowStab;
}

template <>
BlockVector<VectorEp>
StokesAssembler<VectorEp,MatrixEp>::
getZeroVector() const
{
    SHP(VECTOREPETRA) uComp(new VECTOREPETRA(M_velocityFESpace->map(),
                                             LifeV::Unique));

    uComp->zero();

    SHP(VECTOREPETRA) pComp(new VECTOREPETRA(M_pressureFESpace->map(),
                                             LifeV::Unique));

    pComp->zero();

    BlockVector<VectorEp> retVec;
    retVec.resize(M_nComponents);
    retVec.block(0).data() = uComp;
    retVec.block(1).data() = pComp;
    return retVec;
}

template <>
SHP(MATRIXEPETRA)
StokesAssembler<VectorEp, MatrixEp>::
assembleFlowRateJacobian(const GeometricFace& face)
{
    using namespace LifeV;
    using namespace ExpressionAssembly;

    const double dropTolerance(2.0 * std::numeric_limits<double>::min());

    SHP(MAPEPETRA) rangeMap = M_flowRateVectors[face.M_flag]->mapPtr();
    EPETRACOMM comm = rangeMap->commPtr();


    Epetra_Map epetraMap = M_flowRateVectors[face.M_flag]->epetraMap();
    unsigned int numElements = epetraMap.NumMyElements();
    unsigned int numGlobalElements = epetraMap.NumGlobalElements();

    // this should be optimized
    SHP(MATRIXEPETRA) flowRateJacobian;
    flowRateJacobian.reset(new MATRIXEPETRA(M_velocityFESpace->map(), numGlobalElements, false));

    // compute outer product of flowrate vector with itself
    for (unsigned int j = 0; j < numGlobalElements; j++)
    {
        double myvaluecol = 0;

        if (M_flowRateVectors[face.M_flag]->isGlobalIDPresent(j))
            myvaluecol = M_flowRateVectors[face.M_flag]->operator[](j);

        double valuecol = 0;
        comm->SumAll(&myvaluecol, &valuecol, 1);

        if (std::abs(valuecol) > dropTolerance)
        {
            for (unsigned int i = 0; i < numElements; i++)
            {
                unsigned int gdof = epetraMap.GID(i);
                if (M_flowRateVectors[face.M_flag]->isGlobalIDPresent(gdof))
                {
                    double valuerow = M_flowRateVectors[face.M_flag]->operator[](gdof);
                    if (std::abs(valuerow * valuecol) > dropTolerance)
                    {
                        flowRateJacobian->addToCoefficient(gdof, j, valuerow * valuecol);
                    }
                }
            }
        }

    }

    comm->Barrier();

    flowRateJacobian->globalAssemble();

    return flowRateJacobian;
}

template <>
void
StokesAssembler<VectorEp, MatrixEp>::
assembleFlowRateJacobians()
{
    // assemble inflow flow rate vector
    if (M_treeNode->isInletNode())
    {
        auto face = M_treeNode->M_block->getInlet();
        M_flowRateJacobians[face.M_flag].resize(this->M_nComponents,
                                                this->M_nComponents);
        M_flowRateJacobians[face.M_flag].block(0,0).data() = assembleFlowRateJacobian(face);

        apply0DirichletBCsMatrix(M_flowRateJacobians[face.M_flag], 0.0);
    }

    if (M_treeNode->isOutletNode())
    {
        auto faces = M_treeNode->M_block->getOutlets();

        for (auto face : faces)
        {
            M_flowRateJacobians[face.M_flag].resize(this->M_nComponents,
                                                    this->M_nComponents);
            M_flowRateJacobians[face.M_flag].block(0,0).data() = assembleFlowRateJacobian(face);

            apply0DirichletBCsMatrix(M_flowRateJacobians[face.M_flag], 0.0);
        }
    }
}

template <>
BlockVector<FEVECTOR>
StokesAssembler<FEVECTOR, FEMATRIX>::
getLifting(const double& time) const
{
    return getFELifting(time);
}

template <>
void
StokesAssembler<VectorEp, MatrixEp>::
exportSolution(const double& t, const BlockVector<VectorEp>& sol)
{
    *M_velocityExporter = *sol.block(0).data();
    *M_pressureExporter = *sol.block(1).data();

    BlockVector<VectorEp> solCopy(2);
    solCopy.block(0).data() = M_velocityExporter;
    computeFlowRates(solCopy, true);

    exportNorms(t);

    CoutRedirecter ct;
    ct.redirect();
    M_exporter->postProcess(t);
    printlog(CYAN, ct.restore());
}

template <>
MatrixEp
StokesAssembler<VectorEp, MatrixEp>::
getNorm(const unsigned int& fieldIndex)
{
    using namespace LifeV;
    using namespace ExpressionAssembly;

    MatrixEp retMat;
    if (fieldIndex == 0)
    {
        if (!M_massVelocity.data())
        {
            SHP(MATRIXEPETRA) Nu(new MATRIXEPETRA(M_velocityFESpace->map()));

            integrate(elements(M_velocityFESpaceETA->mesh()),
                      M_velocityFESpace->qr(),
                      M_velocityFESpaceETA,
                      M_velocityFESpaceETA,
                      dot(phi_i, phi_j) +
                      dot(grad(phi_i),grad(phi_j))
                  ) >> Nu;

            Nu->globalAssemble();

            BlockMatrix<MatrixEp> normWrap;
            normWrap.resize(1,1);
            normWrap.block(0,0).data() = Nu;

            // note. Applying bcs does not change the norm if Dirichlet bcs are
            // homogeneous (=> lifting) or imposed weakly. Here we impose bcs
            // in order to have the correct conditions in the computation of the
            // supremizers (we have to solve a linear system..)
            apply0DirichletBCsMatrix(normWrap, 1.0);

            M_massVelocity.data() = Nu;
            retMat.data() = Nu;
        }
        else
            retMat = M_massVelocity;
    }
    else
    {
        if (!M_massPressure.data())
        {
            SHP(MATRIXEPETRA) Np(new MATRIXEPETRA(M_pressureFESpace->map()));

            integrate(elements(M_pressureFESpaceETA->mesh()),
                      M_pressureFESpace->qr(),
                      M_pressureFESpaceETA,
                      M_pressureFESpaceETA,
                      phi_i * phi_j
                  ) >> Np;

            Np->globalAssemble();

            M_massPressure.data() = Np;
            retMat.data() = Np;
        }
        else
            retMat = M_massPressure;
    }

    return retMat;
}

template <>
void
StokesAssembler<VectorEp, MatrixEp>::
restrictRBMatrices()
{
    throw new Exception("restrictRBMatrices not available for fem assembler");
}

}
