#include "StokesAssembler.hpp"

namespace RedMA
{

template <>
void
StokesAssembler<VectorEp,MatrixEp>::
assembleStiffness()
{
    using namespace LifeV;
    using namespace LifeV::ExpressionAssembly;

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
}

template <>
void
StokesAssembler<VectorEp,MatrixEp>::
assembleMass()
{

}

}
