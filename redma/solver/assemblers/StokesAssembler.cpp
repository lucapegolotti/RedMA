#include "StokesAssembler.hpp"

namespace RedMA
{

template <>
void
StokesAssembler<VectorEp,MatrixEp>::
assembleStiffness()
{
    using namespace LifeV;
    using namespace ExpressionAssembly;

    M_stiffness.resize(this->M_nComponents, this->M_nComponents);
    bool useFullStrain = M_data("fluid/use_strain", true);

    SHP(MatrixEpetra<double>) A(new MatrixEpetra<double>(M_velocityFESpace->map()));

    if (useFullStrain)
    {
        integrate(elements(M_velocityFESpaceETA->mesh()),
                  M_velocityFESpace->qr(),
                  M_velocityFESpaceETA,
                  M_velocityFESpaceETA,
                  value(0.5 * M_viscosity) *
                  dot(grad(phi_i) + transpose(grad(phi_i)),
                  grad(phi_j) + transpose(grad(phi_j)))
              ) >> A;
    }
    else
    {
        integrate(elements(M_velocityFESpaceETA->mesh()),
                  M_velocityFESpace->qr(),
                  M_velocityFESpaceETA,
                  M_velocityFESpaceETA,
                  value(M_viscosity) *
                  dot(grad(phi_i),grad(phi_j))
              ) >> A;
    }
    A->globalAssemble();

    M_stiffness.block(0,0).data() = A;

    M_bcManager->apply0DirichletMatrix(M_stiffness, getFESpaceBCs(),
                                       getComponentBCs(), 0.0);
}

template <>
void
StokesAssembler<VectorEp,MatrixEp>::
assembleMass()
{
    using namespace LifeV;
    using namespace ExpressionAssembly;

    M_mass.resize(this->M_nComponents, this->M_nComponents);

    SHP(MatrixEpetra<double>) M(new MatrixEpetra<double>(M_velocityFESpace->map()));

    integrate(elements(M_velocityFESpaceETA->mesh()),
              M_velocityFESpace->qr(),
              M_velocityFESpaceETA,
              M_velocityFESpaceETA,
              value(M_density) * dot(phi_i, phi_j)
          ) >> M;

    M->globalAssemble();

    M_mass.block(0,0).data() = M;

    M_bcManager->apply0DirichletMatrix(M_mass, getFESpaceBCs(),
                                       getComponentBCs(), 1.0);
}

template <>
void
StokesAssembler<VectorEp,MatrixEp>::
assembleDivergence()
{
    using namespace LifeV;
    using namespace ExpressionAssembly;

    M_divergence.resize(this->M_nComponents, this->M_nComponents);

    SHP(MatrixEpetra<double>) BT(new MatrixEpetra<double>(M_velocityFESpace->map()));

    integrate(elements(M_velocityFESpaceETA->mesh()),
              M_velocityFESpace->qr(),
              M_velocityFESpaceETA,
              M_pressureFESpaceETA,
              value(-1.0) * phi_j * div(phi_i)
          ) >> BT;

    BT->globalAssemble(M_pressureFESpace->mapPtr(),
                       M_velocityFESpace->mapPtr());

    SHP(MatrixEpetra<double>) B(new MatrixEpetra<double>(M_pressureFESpace->map()));

    integrate(elements(M_velocityFESpaceETA->mesh()),
             M_pressureFESpace->qr(),
             M_pressureFESpaceETA,
             M_velocityFESpaceETA,
             value(-1.0) * phi_i * div(phi_j)
         ) >> B;

    B->globalAssemble(M_velocityFESpace->mapPtr(),
                      M_pressureFESpace->mapPtr());

    M_divergence.block(0,1).data() = BT;
    M_divergence.block(1,0).data() = B;
    M_bcManager->apply0DirichletMatrix(M_divergence, getFESpaceBCs(),
                                       getComponentBCs(), 0.0);
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
SHP(LifeV::VectorEpetra)
StokesAssembler<VectorEp, MatrixEp>::
assembleFlowRateVector(const GeometricFace& face)
{
    using namespace LifeV;
    using namespace ExpressionAssembly;

    SHP(VECTOREPETRA) flowRateVectorRepeated;
    flowRateVectorRepeated.reset(new VECTOREPETRA(M_velocityFESpace->map(),
                                                  Repeated));

    QuadratureBoundary myBDQR(buildTetraBDQR(quadRuleTria7pt));

    integrate(boundary(M_velocityFESpaceETA->mesh(), face.M_flag),
              myBDQR,
              M_velocityFESpaceETA,
              dot(phi_i, Nface)
          ) >> flowRateVectorRepeated;

    flowRateVectorRepeated->globalAssemble();

    SHP(VECTOREPETRA) flowRateVector(new VECTOREPETRA(*flowRateVectorRepeated,
                                                      Unique));
    return flowRateVector;
}

template <>
void
StokesAssembler<VectorEp, MatrixEp>::
assembleFlowRateVectors()
{
    // assemble inflow flow rate vector
    if (M_treeNode->isInletNode())
    {
        auto face = M_treeNode->M_block->getInlet();
        std::cout << "\nInlet normal\n" << std::endl;
        face.print();
        M_flowRateVectors[face.M_flag] = assembleFlowRateVector(face);
    }

    if (M_treeNode->isOutletNode())
    {
        auto faces = M_treeNode->M_block->getOutlets();

        for (auto face : faces)
            M_flowRateVectors[face.M_flag] = assembleFlowRateVector(face);
    }
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
        M_bcManager->apply0DirichletMatrix(M_flowRateJacobians[face.M_flag], getFESpaceBCs(),
                                           getComponentBCs(), 0.0);
    }

    if (M_treeNode->isOutletNode())
    {
        auto faces = M_treeNode->M_block->getOutlets();

        for (auto face : faces)
        {
            M_flowRateJacobians[face.M_flag].resize(this->M_nComponents,
                                                    this->M_nComponents);
            M_flowRateJacobians[face.M_flag].block(0,0).data() = assembleFlowRateJacobian(face);
            M_bcManager->apply0DirichletMatrix(M_flowRateJacobians[face.M_flag], getFESpaceBCs(),
                                               getComponentBCs(), 0.0);
        }
    }
}

template <>
BlockVector<FEVECTOR>
StokesAssembler<FEVECTOR, FEMATRIX>::
getLifting(const double& time) const
{
    BlockVector<FEVECTOR> lifting;
    lifting.resize(2);
    lifting.block(0).data().reset(new VECTOREPETRA(M_velocityFESpace->map(),
                                                   LifeV::Unique));
    lifting.block(0).data()->zero();
    lifting.block(1).data().reset(new VECTOREPETRA(M_pressureFESpace->map(),
                                                   LifeV::Unique));
    lifting.block(1).data()->zero();

    this->M_bcManager->applyDirichletBCs(time, lifting, getFESpaceBCs(),
                                         getComponentBCs());
    return lifting;
}

}
