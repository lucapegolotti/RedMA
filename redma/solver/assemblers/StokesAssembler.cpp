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

    M_stiffness.resize(M_nComponents,M_nComponents);
    bool useFullStrain = M_datafile("fluid/use_strain", true);

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
    M_stiffness.finalize();
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

    M_mass.resize(M_nComponents, M_nComponents);

    SHP(MatrixEpetra<double>) M(new MatrixEpetra<double>(M_velocityFESpace->map()));

    integrate(elements(M_velocityFESpaceETA->mesh()),
              M_velocityFESpace->qr(),
              M_velocityFESpaceETA,
              M_velocityFESpaceETA,
              value(M_density) * dot(phi_i, phi_j)
          ) >> M;

    M->globalAssemble();

    M_mass.block(0,0).data() = M;
    M_mass.finalize();
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

    M_divergence.resize(M_nComponents, M_nComponents);

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
    M_divergence.finalize();
    M_bcManager->apply0DirichletMatrix(M_divergence, getFESpaceBCs(),
                                       getComponentBCs(), 0.0);
}

}
