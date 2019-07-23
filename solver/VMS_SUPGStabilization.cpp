// implementation inspired by StabilizationSUPG.cpp in LifeV
#include "VMS_SUPGStabilization.hpp"

// MACRO TO DEFINE TAU_M
#define TAU_M 	       value(1)/(eval(squareroot,TAU_M_DEN))
#define TAU_M_DEN      (TAU_M_DEN_DT + TAU_M_DEN_VEL + TAU_M_DEN_VISC)
#define TAU_M_DEN_DT   (value(M_density*M_density)*value(M_order*M_order)/value(dt * dt))
#define TAU_M_DEN_VEL  (value(M_density*M_density)*dot(value(M_velocityFESpaceETA, *velocity), G*value(M_velocityFESpaceETA, *velocity)))
#define TAU_M_DEN_VISC (value(M_C_I)*value(M_viscosity*M_viscosity)*dot(G, G))

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
        M_blockMass.reset(new Matrix(M_velocityFESpace->map()));
        M_block00.reset(new Matrix(M_velocityFESpace->map()));
        M_block01.reset(new Matrix(M_velocityFESpace->map()));
        M_block10.reset(new Matrix(M_pressureFESpace->map()));
        M_block11.reset(new Matrix(M_pressureFESpace->map()));
    }

    M_blockMass->zero();
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

    // integrate(elements(M_uFESpace->mesh()),
    //           M_uFESpace->qr(),
    //           M_fespaceUETA,
    //           M_fespaceUETA,
    //           +TAU_M*value(M_density*M_density)*value(M_alpha/M_timestep) * dot( value(M_fespaceUETA, velocity)*grad(phi_i), phi_j)
    //           +TAU_M*value(M_density*M_density)*value(M_alpha/M_timestep) * dot( phi_j*grad(phi_i), value(M_fespaceUETA, velocity))
    //           -TAU_M*value(M_density*M_density)*dot(phi_j*grad(phi_i), value(M_fespaceUETA, velocityRhs))
    //           +TAU_M*value(M_density*M_density)*dot(value(M_fespaceUETA, velocity)*grad(phi_i), phi_j*grad(M_fespaceUETA, velocity))
    //           +TAU_M*value(M_density*M_density)*dot(value(M_fespaceUETA, velocity)*grad(phi_i), value(M_fespaceUETA, velocity)*grad(phi_j))
    //           +TAU_M*value(M_density*M_density)*dot(phi_j*grad(phi_i), value(M_fespaceUETA, velocity)*grad(M_fespaceUETA, velocity))
    //           +TAU_M*value(M_density)*dot( phi_j*grad(phi_i), grad(M_fespacePETA, pressure))
    //           -TAU_M*value(M_density*M_viscosity)*dot(value(M_fespaceUETA, velocity)*grad(phi_i), laplacian(phi_j))
    //           -TAU_M*value(M_density*M_viscosity)*dot(phi_j*grad(phi_i), laplacian(M_fespaceUETA, velocity))
    //           +TAU_C*div(phi_i)*div(phi_j)
    //          ) >> M_block_00;

    integrate(elements(M_velocityFESpace->mesh()),
              M_velocityFESpace->qr(),
              M_velocityFESpaceETA,
              M_velocityFESpaceETA,
               TAU_M*value(M_density*M_density)*dot(value(M_velocityFESpaceETA, *velocity)*grad(phi_i), phi_j)
              +TAU_M*value(M_density*M_density)*dot(phi_j*grad(phi_i), value(M_velocityFESpaceETA, *velocity))
              ) >> M_blockMass;
    M_blockMass->globalAssemble();
}

}  // namespace RedMA
