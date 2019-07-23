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

    if (M_block00 == nullptr && M_block01 == nullptr &&
        M_block10 == nullptr && M_block11 == nullptr)
    {
        M_blockMass00.reset(new Matrix(M_velocityFESpace->map()));
        M_blockMass01.reset(new Matrix(M_velocityFESpace->map()));
        M_blockMass10.reset(new Matrix(M_pressureFESpace->map()));
        M_blockMass11.reset(new Matrix(M_pressureFESpace->map()));
        M_block00.reset(new Matrix(M_velocityFESpace->map()));
        M_block01.reset(new Matrix(M_velocityFESpace->map()));
        M_block10.reset(new Matrix(M_pressureFESpace->map()));
        M_block11.reset(new Matrix(M_pressureFESpace->map()));
    }

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
               TAU_M*value(M_density*M_density)*dot(value(M_velocityFESpaceETA, *velocity)*grad(phi_i), phi_j)
              +TAU_M*value(M_density*M_density)*dot(phi_j*grad(phi_i), value(M_velocityFESpaceETA, *velocity))
          ) >> M_blockMass00;
    M_blockMass00->globalAssemble();

    integrate(elements(M_velocityFESpace->mesh()),
              M_velocityFESpace->qr(),
              M_velocityFESpaceETA,
              M_velocityFESpaceETA,
               TAU_M*value(M_density*M_density)*dot(value(M_velocityFESpaceETA, *velocity)*grad(phi_i), phi_j*grad(M_velocityFESpaceETA, *velocity))
              +TAU_M*value(M_density*M_density)*dot(value(M_velocityFESpaceETA, *velocity)*grad(phi_i), value(M_velocityFESpaceETA, *velocity)*grad(phi_j))
              -TAU_M*value(M_density*M_viscosity)*dot(value(M_velocityFESpaceETA, *velocity)*grad(phi_i), laplacian(phi_j))
              -TAU_M*value(M_density*M_density)*dot(phi_j*grad(phi_i), value(M_velocityFESpaceETA, *velocityRhs))
              +TAU_M*value(M_density*M_density)*dot(phi_j*grad(phi_i), value(M_velocityFESpaceETA, *velocity)*grad(M_velocityFESpaceETA, *velocity))
              +TAU_M*value(M_density)*dot(phi_j*grad(phi_i), grad(M_pressureFESpaceETA, *pressure))
              -TAU_M*value(M_density*M_viscosity)*dot(phi_j*grad(phi_i), laplacian(M_velocityFESpaceETA, *velocity))
              +TAU_C*div(phi_i)*div(phi_j)
              ) >> M_block00;
    M_block00->globalAssemble();

    integrate(elements(M_velocityFESpace->mesh()),
              M_pressureFESpace->qr(),
              M_pressureFESpaceETA,
              M_velocityFESpaceETA,
              TAU_M*value(M_density)*dot( grad(phi_i), phi_j)
             ) >> M_blockMass10;
    M_blockMass10->globalAssemble(M_pressureFESpace->mapPtr(),
                                  M_velocityFESpace->mapPtr());

    integrate(elements(M_velocityFESpace->mesh()),
              M_pressureFESpace->qr(),
              M_pressureFESpaceETA,
              M_velocityFESpaceETA,
               TAU_M*value(M_density)*dot( grad(phi_i), phi_j*grad(M_velocityFESpaceETA, *velocity))
              +TAU_M*value(M_density)*dot( grad(phi_i), value(M_velocityFESpaceETA, *velocity)*grad(phi_j))
              -TAU_M*value(M_viscosity)*dot(grad(phi_i), laplacian(phi_j))
              ) >> M_block10;
    M_block10->globalAssemble(M_pressureFESpace->mapPtr(),
                              M_velocityFESpace->mapPtr());

    integrate(elements(M_velocityFESpace->mesh()),
              M_velocityFESpace->qr(),
              M_velocityFESpaceETA,
              M_pressureFESpaceETA,
              TAU_M*value(M_density)*dot(value(M_velocityFESpaceETA, *velocity)*grad(phi_i), grad(phi_j))
              ) >> M_block01;
    M_block01->globalAssemble(M_pressureFESpace->mapPtr(),
                              M_velocityFESpace->mapPtr());

    integrate(elements(M_velocityFESpace->mesh()),
              M_pressureFESpace->qr(),
              M_pressureFESpaceETA,
              M_pressureFESpaceETA,
              TAU_M*dot(grad(phi_i), grad(phi_j))
              ) >> M_block11;
    M_block11->globalAssemble();
}

}  // namespace RedMA
