// implementation inspired by StabilizationSUPG.cpp in LifeV
#include "VMS_SUPGStabilization.hpp"

// MACRO TO DEFINE TAU_M
#define TAU_M 	       value(1)/(eval(squareroot,TAU_M_DEN))
#define TAU_M_DEN      (TAU_M_DEN_DT + TAU_M_DEN_VEL + TAU_M_DEN_VISC)
#define TAU_M_DEN_DT   (value(M_density*M_density)*value(M_order*M_order)/value(dt * dt))
#define TAU_M_DEN_VEL  (value(M_density*M_density)*dot(value(M_velocityFESpaceETA, *velocity), G*value(M_velocityFESpaceETA, *velocity)))
#define TAU_M_DEN_VISC (value(M_C_I)*value(M_viscosity*M_viscosity)*dot(G, G))

#define TAU_C ( value(1.0)/( dot(g, TAU_M*g ) ) )

namespace RedMA
{

VMS_SUPGStabilization::
VMS_SUPGStabilization(unsigned int         order,
                     FESpacePtr           fespaceVelocity,
                     FESpacePtr           fespacePressure,
                     ETFESpaceVelocityPtr etfespaceVelocity,
                     ETFESpacePressurePtr etfespacePressure) :
  M_order(order),
  M_velocityFESpace(fespaceVelocity),
  M_pressureFESpace(fespacePressure),
  M_velocityFESpaceETA(etfespaceVelocity),
  M_pressureFESpaceETA(etfespacePressure)
{
    if (order == 1)
        M_C_I = 30;
    else if (order == 2)
        M_C_I = 60;
    else if (order == 3)
        M_C_I = 120;
    else if (order == 4)
        M_C_I = 240;
    else
    {
        std::string msg = "Please implement a suitable value for ";
        msg += " M_C_I for your velocity FE order";
        throw Exception(msg);
    }
}

void
VMS_SUPGStabilization::
setDensityAndViscosity(double density, double viscosity)
{
    M_density = density;
    M_viscosity = viscosity;
}

void
VMS_SUPGStabilization::
assembleBlocks(VectorPtr velocity, VectorPtr pressure,
               VectorPtr velocityRhs, double dt)
{
    using namespace LifeV::ExpressionAssembly;

    if (M_block00Jac == nullptr && M_block01Jac == nullptr &&
        M_block10Jac == nullptr && M_block11Jac == nullptr)
    {
        M_blockMass00Jac.reset(new Matrix(M_velocityFESpace->map()));
        M_blockMass01Jac.reset(new Matrix(M_velocityFESpace->map()));
        M_blockMass10Jac.reset(new Matrix(M_pressureFESpace->map()));
        M_blockMass11Jac.reset(new Matrix(M_pressureFESpace->map()));
        M_block00Jac.reset(new Matrix(M_velocityFESpace->map()));
        M_block01Jac.reset(new Matrix(M_velocityFESpace->map()));
        M_block10Jac.reset(new Matrix(M_pressureFESpace->map()));
        M_block11Jac.reset(new Matrix(M_pressureFESpace->map()));
        M_blockMass00.reset(new Matrix(M_velocityFESpace->map()));
        M_blockMass01.reset(new Matrix(M_velocityFESpace->map()));
        M_blockMass10.reset(new Matrix(M_pressureFESpace->map()));
        M_blockMass11.reset(new Matrix(M_pressureFESpace->map()));
        M_block00.reset(new Matrix(M_velocityFESpace->map()));
        M_block01.reset(new Matrix(M_velocityFESpace->map()));
        M_block10.reset(new Matrix(M_pressureFESpace->map()));
        M_block11.reset(new Matrix(M_pressureFESpace->map()));
    }

    M_blockMass00Jac->zero();
    M_blockMass01Jac->zero();
    M_blockMass10Jac->zero();
    M_blockMass11Jac->zero();
    M_block00Jac->zero();
    M_block01Jac->zero();
    M_block10Jac->zero();
    M_block11Jac->zero();
    M_blockMass00->zero();
    M_blockMass01->zero();
    M_blockMass10->zero();
    M_blockMass11->zero();
    M_block00->zero();
    M_block01->zero();
    M_block10->zero();
    M_block11->zero();

    std::shared_ptr<LifeV::SquareRoot> squareroot(new LifeV::SquareRoot());

    LifeV::MatrixSmall<3, 3> Eye;
    Eye *= 0.0;
    Eye[0][0] = 1;
    Eye[1][1] = 1;
    Eye[2][2] = 1;

    integrate(elements(M_velocityFESpace->mesh()),
              M_velocityFESpace->qr(),
              M_velocityFESpaceETA,
              M_velocityFESpaceETA,
               TAU_M * value(M_density*M_density) *
                       dot(value(M_velocityFESpaceETA, *velocity) *
                       grad(phi_i), phi_j)
              +TAU_M * value(M_density*M_density) *
                       dot(phi_j*grad(phi_i),
                       value(M_velocityFESpaceETA, *velocity))
          ) >> M_blockMass00Jac;

    integrate(elements(M_velocityFESpace->mesh()),
              M_velocityFESpace->qr(),
              M_velocityFESpaceETA,
              M_velocityFESpaceETA,
               TAU_M * value(M_density*M_density) *
                       dot(value(M_velocityFESpaceETA, *velocity) *
                       grad(phi_i), phi_j)
          ) >> M_blockMass00;

    integrate(elements(M_velocityFESpace->mesh()),
              M_velocityFESpace->qr(),
              M_velocityFESpaceETA,
              M_velocityFESpaceETA,
               TAU_M * value(M_density*M_density) *
                       dot(value(M_velocityFESpaceETA, *velocity) *
                       grad(phi_i), phi_j*grad(M_velocityFESpaceETA, *velocity))
              +TAU_M * value(M_density*M_density) *
                       dot(value(M_velocityFESpaceETA, *velocity) *
                       grad(phi_i), value(M_velocityFESpaceETA, *velocity) *
                       grad(phi_j))
              -TAU_M * value(M_density*M_viscosity) *
                       dot(value(M_velocityFESpaceETA, *velocity) *
                       grad(phi_i), laplacian(phi_j))
              -TAU_M * value(M_density*M_density) *
                       dot(phi_j*grad(phi_i), value(M_velocityFESpaceETA, *velocityRhs))
              +TAU_M * value(M_density*M_density) *
                       dot(phi_j*grad(phi_i), value(M_velocityFESpaceETA, *velocity) *
                       grad(M_velocityFESpaceETA, *velocity))
              +TAU_M * value(M_density) *
                       dot(phi_j*grad(phi_i), grad(M_pressureFESpaceETA, *pressure))
              -TAU_M * value(M_density*M_viscosity) *
                       dot(phi_j*grad(phi_i), laplacian(M_velocityFESpaceETA, *velocity))
              +TAU_C * div(phi_i)*div(phi_j)
              ) >> M_block00Jac;

    integrate(elements(M_velocityFESpace->mesh()),
              M_velocityFESpace->qr(),
              M_velocityFESpaceETA,
              M_velocityFESpaceETA,
               TAU_M * value(M_density*M_density) *
                       dot(value(M_velocityFESpaceETA, *velocity) *
                       grad(phi_i), phi_j*grad(M_velocityFESpaceETA, *velocity))
              +TAU_M * value(M_density*M_density) *
                       dot(value(M_velocityFESpaceETA, *velocity) *
                       grad(phi_i), value(M_velocityFESpaceETA, *velocity) *
                       grad(phi_j))
              -TAU_M * value(M_density*M_viscosity) *
                       dot(value(M_velocityFESpaceETA, *velocity) *
                       grad(phi_i), laplacian(phi_j))
              -TAU_M * value(M_density*M_density) *
                       dot(phi_j*grad(phi_i), value(M_velocityFESpaceETA, *velocityRhs))
              +TAU_C * div(phi_i)*div(phi_j)
              ) >> M_block00Jac;

    integrate(elements(M_velocityFESpace->mesh()),
              M_pressureFESpace->qr(),
              M_pressureFESpaceETA,
              M_velocityFESpaceETA,
              TAU_M * value(M_density)*dot(grad(phi_i), phi_j)
              ) >> M_blockMass10Jac;

    integrate(elements(M_velocityFESpace->mesh()),
           M_pressureFESpace->qr(),
           M_pressureFESpaceETA,
           M_velocityFESpaceETA,
           TAU_M * value(M_density)*dot(grad(phi_i), phi_j)
              ) >> M_blockMass10;

    integrate(elements(M_velocityFESpace->mesh()),
              M_pressureFESpace->qr(),
              M_pressureFESpaceETA,
              M_velocityFESpaceETA,
               TAU_M * value(M_density) *
                       dot(grad(phi_i), phi_j *
                       grad(M_velocityFESpaceETA, *velocity))
              +TAU_M * value(M_density) *
                       dot(grad(phi_i), value(M_velocityFESpaceETA, *velocity) *
                       grad(phi_j))
              -TAU_M * value(M_viscosity) * dot(grad(phi_i), laplacian(phi_j))
              ) >> M_block10Jac;

    integrate(elements(M_velocityFESpace->mesh()),
              M_pressureFESpace->qr(),
              M_pressureFESpaceETA,
              M_velocityFESpaceETA,
              TAU_M * value(M_density) *
                      dot(grad(phi_i), phi_j *
                      grad(M_velocityFESpaceETA, *velocity))
             -TAU_M * value(M_viscosity) * dot(grad(phi_i), laplacian(phi_j))
              ) >> M_block10;

    integrate(elements(M_velocityFESpace->mesh()),
              M_velocityFESpace->qr(),
              M_velocityFESpaceETA,
              M_pressureFESpaceETA,
              TAU_M * value(M_density) *
                      dot(value(M_velocityFESpaceETA, *velocity) *
                      grad(phi_i), grad(phi_j))
              ) >> M_block01Jac;

    integrate(elements(M_velocityFESpace->mesh()),
           M_velocityFESpace->qr(),
           M_velocityFESpaceETA,
           M_pressureFESpaceETA,
           TAU_M * value(M_density) *
                   dot(value(M_velocityFESpaceETA, *velocity) *
                   grad(phi_i), grad(phi_j))
           ) >> M_block01;

    integrate(elements(M_velocityFESpace->mesh()),
              M_pressureFESpace->qr(),
              M_pressureFESpaceETA,
              M_pressureFESpaceETA,
              TAU_M * dot(grad(phi_i), grad(phi_j))
              ) >> M_block11Jac;

    M_blockMass00Jac->globalAssemble();
    M_blockMass00->globalAssemble();
    M_block00Jac->globalAssemble();
    M_block00->globalAssemble();
    M_blockMass10Jac->globalAssemble(M_velocityFESpace->mapPtr(),
                                     M_pressureFESpace->mapPtr());
    M_blockMass10->globalAssemble(M_velocityFESpace->mapPtr(),
                                  M_pressureFESpace->mapPtr());
    M_block10Jac->globalAssemble(M_velocityFESpace->mapPtr(),
                               M_pressureFESpace->mapPtr());
    M_block10->globalAssemble(M_velocityFESpace->mapPtr(),
                              M_pressureFESpace->mapPtr());
    M_block01Jac->globalAssemble(M_pressureFESpace->mapPtr(),
                                 M_velocityFESpace->mapPtr());
    M_block01->globalAssemble(M_pressureFESpace->mapPtr(),
                              M_velocityFESpace->mapPtr());
    M_blockMass01Jac->globalAssemble(M_pressureFESpace->mapPtr(),
                                     M_velocityFESpace->mapPtr());
    M_blockMass01->globalAssemble(M_pressureFESpace->mapPtr(),
                                  M_velocityFESpace->mapPtr());
    M_block11Jac->globalAssemble();
    M_block11->globalAssemble();
    M_blockMass11Jac->globalAssemble();
    M_blockMass11->globalAssemble();

    // M_blockMass00Jac->zero();
    // M_blockMass01Jac->zero();
    // M_blockMass10Jac->zero();
    // M_blockMass11Jac->zero();
    // M_block00Jac->zero();
    // M_block01Jac->zero();
    // M_block10Jac->zero();
    // M_block11Jac->zero();
    // M_blockMass00->zero();
    // M_blockMass01->zero();
    // M_blockMass10->zero();
    // M_blockMass11->zero();
    // M_block00->zero();
    // M_block01->zero();
    // M_block10->zero();
    // M_block11->zero();
}

}  // namespace RedMA
