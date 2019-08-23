// implementation inspired by StabilizationSUPG.cpp in LifeV
#include "VMS_SUPGStabilization.hpp"

// MACRO TO DEFINE TAU_M
// #define TAU_M 	       value(1.0)/(eval(squareroot,TAU_M_DEN))
// #define TAU_M_DEN      (TAU_M_DEN_DT + TAU_M_DEN_VEL + TAU_M_DEN_VISC)
// #define TAU_M_DEN_DT   (value(M_density*M_density)*value(M_timeOrder*M_timeOrder)/value(dt * dt))
// #define TAU_M_DEN_VEL  (value(M_density*M_density)*dot(value(M_velocityFESpaceETA, *velocityRep), G*value(M_velocityFESpaceETA, *velocityRep)))
// #define TAU_M_DEN_VISC (value(M_C_I)*value(M_viscosity*M_viscosity)*dot(G, G))
//
// #define TAU_C ( value(1.0)/( dot(g, TAU_M*g ) ) )

#define TAU_M 	       value(1.0)/(eval(squareroot,TAU_M_DEN))
#define TAU_M_DEN      (TAU_M_DEN_DT + TAU_M_DEN_VEL + TAU_M_DEN_VISC)
#define TAU_M_DEN_DT   (value(M_density*M_density)*value(M_timeOrder*M_timeOrder)/value(dt*dt))
#define TAU_M_DEN_VEL  (value(M_density*M_density)*dot(value(M_velocityFESpaceETA, *velocityRep)/h_K, value(M_velocityFESpaceETA, *velocityRep)/h_K))
#define TAU_M_DEN_VISC (value(M_C_I)*value(M_viscosity*M_viscosity)/(h_K*h_K*h_K*h_K))

#define TAU_C h_K*h_K/(TAU_M)

namespace RedMA
{

VMS_SUPGStabilization::
VMS_SUPGStabilization(unsigned int         timeOrder,
                      unsigned int         velocityOrder,
                      FESpacePtr           fespaceVelocity,
                      FESpacePtr           fespacePressure,
                      ETFESpaceVelocityPtr etfespaceVelocity,
                      ETFESpacePressurePtr etfespacePressure) :
  M_timeOrder(timeOrder),
  M_velocityOrder(velocityOrder),
  M_velocityFESpace(fespaceVelocity),
  M_pressureFESpace(fespacePressure),
  M_velocityFESpaceETA(etfespaceVelocity),
  M_pressureFESpaceETA(etfespacePressure)
{
    if (velocityOrder == 1)
        M_C_I = 30;
    else if (velocityOrder == 2)
        M_C_I = 60;
    else if (velocityOrder == 3)
        M_C_I = 120;
    else if (velocityOrder == 4)
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

    VectorPtr velocityRep(new Vector(*velocity, LifeV::Repeated));
	VectorPtr velocityRhsRep(new Vector(*velocityRhs, LifeV::Repeated));
	VectorPtr pressureRep(new Vector(*pressure, LifeV::Repeated));

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

    std::shared_ptr<LifeV::SquareRoot> squareroot(new LifeV::SquareRoot());

    integrate(elements(M_velocityFESpace->mesh()),
              M_velocityFESpace->qr(),
              M_velocityFESpaceETA,
              M_velocityFESpaceETA,
              TAU_M * value(M_density*M_density)
                    * dot(value(M_velocityFESpaceETA, *velocityRep)
                    * grad(phi_i), phi_j)
              ) >> M_blockMass00;

    // note that in Davide's thesis, the stabilization term for the pressure
    // has no 1/density in front of the test whereas on wikipedia there is
    // https://en.wikipedia.org/wiki/Streamline_upwind_Petrov%E2%80%93Galerkin_pressure-stabilizing_Petrov%E2%80%93Galerkin_formulation_for_incompressible_Navier%E2%80%93Stokes_equations

    integrate(elements(M_velocityFESpace->mesh()),
              M_pressureFESpace->qr(),
              M_pressureFESpaceETA,
              M_velocityFESpaceETA,
              TAU_M * value(M_density) * dot(grad(phi_i), phi_j)
              ) >> M_blockMass10;

    integrate(elements(M_velocityFESpace->mesh()),
              M_velocityFESpace->qr(),
              M_velocityFESpaceETA,
              M_velocityFESpaceETA,
              TAU_M * (value(M_density*M_density)
                       * dot(value(M_velocityFESpaceETA, *velocityRep)
                       * grad(phi_i), phi_j)
                    +  value(M_density*M_density)
                       * dot(phi_j * grad(phi_i),
                         value(M_velocityFESpaceETA, *velocityRep)))
             ) >> M_blockMass00Jac;

    integrate(elements(M_velocityFESpace->mesh()),
              M_pressureFESpace->qr(),
              M_pressureFESpaceETA,
              M_velocityFESpaceETA,
              TAU_M * value(M_density) * dot(grad(phi_i), phi_j)
              ) >> M_blockMass10Jac;

    integrate(elements(M_velocityFESpace->mesh()),
              M_velocityFESpace->qr(),
              M_velocityFESpaceETA,
              M_velocityFESpaceETA,
              TAU_M * (dot(value(M_velocityFESpaceETA, *velocityRep) * grad(phi_i),
                       value(M_density*M_density) *
                       phi_j * grad(M_velocityFESpaceETA, *velocityRep)
                    +  value(M_density*M_density) *
                       value(M_velocityFESpaceETA, *velocityRep) * grad(phi_j)
                    -  value(M_density*M_viscosity) *
                       laplacian(phi_j))
                    +  dot(phi_j * grad(phi_i),
                       value(M_density) *
                       value(M_velocityFESpaceETA, *velocityRhsRep)
                    +  value(M_density*M_density) *
                       value(M_velocityFESpaceETA, *velocityRep) *
                       grad(M_velocityFESpaceETA, *velocityRep)
                    +  value(M_density) *
                       grad(M_pressureFESpaceETA, *pressureRep)
                    -  value(M_density*M_viscosity) *
                       laplacian(M_velocityFESpaceETA, *velocityRep)))
            + TAU_C *  div(phi_i)*div(phi_j)
              ) >> M_block00Jac;

    integrate(elements(M_velocityFESpace->mesh()),
              M_pressureFESpace->qr(),
              M_pressureFESpaceETA,
              M_velocityFESpaceETA,
              TAU_M * (dot(grad(phi_i),
                       value(M_density) *
                       phi_j * grad(M_velocityFESpaceETA, *velocity)
                    +  value(M_density) *
                       value(M_velocityFESpaceETA, *velocityRep) * grad(phi_j)
                    -  value(M_viscosity) *
                       laplacian(phi_j)))
              ) >> M_block10Jac;

    integrate(elements(M_velocityFESpace->mesh()),
              M_velocityFESpace->qr(),
              M_velocityFESpaceETA,
              M_pressureFESpaceETA,
              TAU_M * (value(M_density)
                       * dot(value(M_velocityFESpaceETA, *velocityRep)
                       * grad(phi_i), grad(phi_j)))
              ) >> M_block01Jac;

    integrate(elements(M_velocityFESpace->mesh()),
              M_pressureFESpace->qr(),
              M_pressureFESpaceETA,
              M_pressureFESpaceETA,
              TAU_M * dot(grad(phi_i), grad(phi_j))
              ) >> M_block11Jac;

    M_blockMass00Jac->globalAssemble();
    M_blockMass00->globalAssemble();
    M_block00Jac->globalAssemble();
    M_blockMass10Jac->globalAssemble(M_velocityFESpace->mapPtr(),
                                     M_pressureFESpace->mapPtr());
    M_blockMass10->globalAssemble(M_velocityFESpace->mapPtr(),
                                  M_pressureFESpace->mapPtr());
    M_block10Jac->globalAssemble(M_velocityFESpace->mapPtr(),
                               M_pressureFESpace->mapPtr());
    M_block01Jac->globalAssemble(M_pressureFESpace->mapPtr(),
                                 M_velocityFESpace->mapPtr());
    M_blockMass01Jac->globalAssemble(M_pressureFESpace->mapPtr(),
                                     M_velocityFESpace->mapPtr());
    M_blockMass01->globalAssemble(M_pressureFESpace->mapPtr(),
                                  M_velocityFESpace->mapPtr());
    M_block11Jac->globalAssemble();
    M_blockMass11Jac->globalAssemble();
    M_blockMass11->globalAssemble();
}

VMS_SUPGStabilization::MatrixPtr
VMS_SUPGStabilization::
assembleMassWithVelocity(VectorPtr velocity, double dt)
{
    using namespace LifeV::ExpressionAssembly;

    VectorPtr velocityRep(new Vector(*velocity, LifeV::Repeated));

    MatrixPtr retMatrix;
    retMatrix.reset(new Matrix(M_velocityFESpace->map()));
    retMatrix->zero();
    std::shared_ptr<LifeV::SquareRoot> squareroot(new LifeV::SquareRoot());

    integrate(elements(M_velocityFESpace->mesh()),
              M_velocityFESpace->qr(),
              M_velocityFESpaceETA,
              M_velocityFESpaceETA,
              TAU_M * value(M_density*M_density)
                    * dot(phi_j*grad(phi_i),
                      value(M_velocityFESpaceETA, *velocityRep))
              ) >> retMatrix;
    retMatrix->globalAssemble();

    return retMatrix;
}

VMS_SUPGStabilization::VectorPtr
VMS_SUPGStabilization::
velocityResidual(VectorPtr velocity, VectorPtr pressure,
                 VectorPtr velocityRhs, double dt)
{
    using namespace LifeV::ExpressionAssembly;

    VectorPtr velocityRep(new Vector(*velocity, LifeV::Repeated));
    VectorPtr velocityRhsRep(new Vector(*velocityRhs, LifeV::Repeated));
    VectorPtr pressureRep(new Vector(*pressure, LifeV::Repeated));

    std::shared_ptr<LifeV::SquareRoot> squareroot(new LifeV::SquareRoot());

    VectorPtr retVec(new Vector(M_velocityFESpace->map()));
    retVec->zero();

    integrate(elements(M_velocityFESpace->mesh()),
              M_velocityFESpace->qr(),
              M_velocityFESpaceETA,
              TAU_M * (dot(value(M_velocityFESpaceETA, *velocityRep) * grad(phi_i),
                       value(M_density*M_density) *
                       value(M_velocityFESpaceETA, *velocityRep) *
                       grad(M_velocityFESpaceETA, *velocityRep)
                    -  value(M_density*M_density) *
                       value(M_velocityFESpaceETA, *velocityRhsRep)
                    +  value(M_density) *
                       grad(M_pressureFESpaceETA, *pressureRep)
                    -  value(M_density*M_viscosity) *
                       laplacian(M_velocityFESpaceETA, *velocityRep)))
              +TAU_C * div(phi_i) * trace(grad(M_velocityFESpaceETA, *velocityRep))
            ) >> retVec;

    return retVec;
}

VMS_SUPGStabilization::VectorPtr
VMS_SUPGStabilization::
pressureResidual(VectorPtr velocity, VectorPtr pressure,
                 VectorPtr velocityRhs, double dt)
{
    using namespace LifeV::ExpressionAssembly;

    VectorPtr velocityRep(new Vector(*velocity, LifeV::Repeated));
    VectorPtr velocityRhsRep(new Vector(*velocityRhs, LifeV::Repeated));
    VectorPtr pressureRep(new Vector(*pressure, LifeV::Repeated));

    std::shared_ptr<LifeV::SquareRoot> squareroot(new LifeV::SquareRoot());

    VectorPtr retVec(new Vector(M_pressureFESpace->map()));
    retVec->zero();

    integrate(elements(M_velocityFESpace->mesh()),
              M_pressureFESpace->qr(),
              M_pressureFESpaceETA,
               TAU_M * (dot(grad(phi_i),
                        value(M_density) *
                        value(M_velocityFESpaceETA, *velocityRep) *
                        grad(M_velocityFESpaceETA, *velocityRep)
                     -  value(M_density) *
                        value(M_velocityFESpaceETA, *velocityRhsRep)
                     +  grad(M_pressureFESpaceETA, *pressureRep)
                     -  value(M_viscosity) *
                        laplacian(M_velocityFESpaceETA, *velocityRep)))
             ) >> retVec;

    return retVec;
}


}  // namespace RedMA
