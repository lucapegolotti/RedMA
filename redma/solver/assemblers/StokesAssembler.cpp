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
BlockVector<VectorEp>
StokesAssembler<VectorEp,MatrixEp>::
getZeroVector() const
{
    BlockVector<VectorEp> retVec;
    retVec.resize(M_nComponents);

    if (isOwned())
    {
        SHP(VECTOREPETRA) uComp(new VECTOREPETRA(M_velocityFESpace->map(),
                                                 LifeV::Unique));

        uComp->zero();

        SHP(VECTOREPETRA) pComp(new VECTOREPETRA(M_pressureFESpace->map(),
                                                 LifeV::Unique));

        pComp->zero();

        retVec.block(0).data() = uComp;
        retVec.block(1).data() = pComp;
    }
    return retVec;
}

}
