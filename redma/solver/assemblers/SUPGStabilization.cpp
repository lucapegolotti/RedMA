// implementation inspired by StabilizationSUPG.cpp in LifeV
#include "SUPGStabilization.hpp"

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
#define TAU_M_DEN_DT   (value(M_density*M_density)*value(M_timeOrder*M_timeOrder)/value(M_dt*M_dt))
#define TAU_M_DEN_VEL  (value(M_density*M_density)*dot(value(M_velocityFESpaceETA, *velocityRep)/h_K, value(M_velocityFESpaceETA, *velocityRep)/h_K))
#define TAU_M_DEN_VISC (value(M_C_I)*value(M_viscosity*M_viscosity)/(h_K*h_K*h_K*h_K))

#define TAU_C h_K*h_K/(TAU_M)

namespace RedMA
{

SUPGStabilization::
SUPGStabilization(const DataContainer& data,
                  SHP(FESPACE) fespaceVelocity,
                  SHP(FESPACE) fespacePressure,
                  SHP(ETFESPACE3) etfespaceVelocity,
                  SHP(ETFESPACE1) etfespacePressure) :
  M_timeOrder(data("time_discretization/order", 2)),
  M_dt(data("time_discretization/dt", 0.01)),
  M_velocityFESpace(fespaceVelocity),
  M_pressureFESpace(fespacePressure),
  M_velocityFESpaceETA(etfespaceVelocity),
  M_pressureFESpaceETA(etfespacePressure)
{
    std::string velocityOrder = data("fluid/velocity_order", "P1");
    if (!std::strcmp(velocityOrder.c_str(),"P1"))
        M_C_I = 30;
    else if (!std::strcmp(velocityOrder.c_str(),"P2"))
        M_C_I = 60;
    else if (!std::strcmp(velocityOrder.c_str(),"P3"))
        M_C_I = 120;
    else if (!std::strcmp(velocityOrder.c_str(),"P4"))
        M_C_I = 240;
    else
    {
        std::string msg = "Please implement a suitable value for ";
        msg += " M_C_I for your velocity FE order";
        throw Exception(msg);
    }
}

void
SUPGStabilization::
setDensityAndViscosity(const double& density, const double& viscosity)
{
    M_density = density;
    M_viscosity = viscosity;
}

BlockMatrix<MatrixEp>
SUPGStabilization::
getMass(const BlockVector<VectorEp>& sol, const BlockVector<VectorEp>& rhs)
{
    using namespace LifeV;
    using namespace LifeV::ExpressionAssembly;
    SHP(SquareRoot) squareroot(new SquareRoot());

    SHP(VECTOREPETRA) velocityRep(new VECTOREPETRA(*sol.block(0).data(), Repeated));
    SHP(VECTOREPETRA) velocityRhsRep(new VECTOREPETRA(*rhs.block(0).data(), Repeated));
    SHP(VECTOREPETRA) pressureRep(new VECTOREPETRA(*sol.block(1).data(), Repeated));

    SHP(MATRIXEPETRA) mass00(new MATRIXEPETRA(M_velocityFESpace->map()));
    SHP(MATRIXEPETRA) mass01(new MATRIXEPETRA(M_velocityFESpace->map()));
    SHP(MATRIXEPETRA) mass10(new MATRIXEPETRA(M_pressureFESpace->map()));
    SHP(MATRIXEPETRA) mass11(new MATRIXEPETRA(M_pressureFESpace->map()));

    mass00->zero();
    mass01->zero();
    mass10->zero();
    mass11->zero();

    integrate(elements(M_velocityFESpace->mesh()),
              M_velocityFESpace->qr(),
              M_velocityFESpaceETA,
              M_velocityFESpaceETA,
              TAU_M * value(M_density*M_density)
                    * dot(value(M_velocityFESpaceETA, *velocityRep)
                    * grad(phi_i), phi_j)
              ) >> mass00;

    // note that in Davide's thesis, the stabilization term for the pressure
    // has no 1/density in front of the test whereas on wikipedia there is
    // https://en.wikipedia.org/wiki/Streamline_upwind_Petrov%E2%80%93Galerkin_pressure-stabilizing_Petrov%E2%80%93Galerkin_formulation_for_incompressible_Navier%E2%80%93Stokes_equations

    integrate(elements(M_velocityFESpace->mesh()),
              M_pressureFESpace->qr(),
              M_pressureFESpaceETA,
              M_velocityFESpaceETA,
              TAU_M * value(M_density) * dot(grad(phi_i), phi_j)
              ) >> mass10;

    mass00->globalAssemble();

    mass01->globalAssemble(M_pressureFESpace->mapPtr(),
                           M_velocityFESpace->mapPtr());

    mass10->globalAssemble(M_velocityFESpace->mapPtr(),
                           M_pressureFESpace->mapPtr());

    mass11->globalAssemble();

    M_mass.resize(2,2);
    M_mass.block(0,0).data() = mass00;
    M_mass.block(0,1).data() = mass01;
    M_mass.block(1,0).data() = mass10;
    M_mass.block(1,1).data() = mass11;

    return M_mass;
}

BlockMatrix<MatrixEp>
SUPGStabilization::
getMassJac(const BlockVector<VectorEp>& sol, const BlockVector<VectorEp>& rhs)
{
    using namespace LifeV;
    using namespace LifeV::ExpressionAssembly;
    SHP(SquareRoot) squareroot(new SquareRoot());

    SHP(VECTOREPETRA) velocityRep(new VECTOREPETRA(*sol.block(0).data(), Repeated));
    SHP(VECTOREPETRA) velocityRhsRep(new VECTOREPETRA(*rhs.block(0).data(), Repeated));
    SHP(VECTOREPETRA) pressureRep(new VECTOREPETRA(*sol.block(1).data(), Repeated));

    SHP(MATRIXEPETRA) massjac00(new MATRIXEPETRA(M_velocityFESpace->map()));
    SHP(MATRIXEPETRA) massjac01(new MATRIXEPETRA(M_velocityFESpace->map()));
    SHP(MATRIXEPETRA) massjac10(new MATRIXEPETRA(M_pressureFESpace->map()));
    SHP(MATRIXEPETRA) massjac11(new MATRIXEPETRA(M_pressureFESpace->map()));

    massjac00->zero();
    massjac01->zero();
    massjac10->zero();
    massjac11->zero();

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
             ) >> massjac00;

    integrate(elements(M_velocityFESpace->mesh()),
              M_pressureFESpace->qr(),
              M_pressureFESpaceETA,
              M_velocityFESpaceETA,
              TAU_M * value(M_density) * dot(grad(phi_i), phi_j)
              ) >> massjac10;

    massjac00->globalAssemble();

    massjac10->globalAssemble(M_velocityFESpace->mapPtr(),
                              M_pressureFESpace->mapPtr());

    massjac01->globalAssemble(M_pressureFESpace->mapPtr(),
                              M_velocityFESpace->mapPtr());

    massjac11->globalAssemble();

    M_massJac.resize(2,2);
    M_massJac.block(0,0).data() = massjac00;
    M_massJac.block(0,1).data() = massjac01;
    M_massJac.block(1,0).data() = massjac10;
    M_massJac.block(1,1).data() = massjac11;

    return M_massJac;
}

BlockMatrix<MatrixEp>
SUPGStabilization::
getJac(const BlockVector<VectorEp>& sol, const BlockVector<VectorEp>& rhs)
{
    using namespace LifeV;
    using namespace LifeV::ExpressionAssembly;
    SHP(SquareRoot) squareroot(new SquareRoot());

    SHP(VECTOREPETRA) velocityRep(new VECTOREPETRA(*sol.block(0).data(), Repeated));
    SHP(VECTOREPETRA) velocityRhsRep(new VECTOREPETRA(*rhs.block(0).data(), Repeated));
    SHP(VECTOREPETRA) pressureRep(new VECTOREPETRA(*sol.block(1).data(), Repeated));

    SHP(MATRIXEPETRA) jac00(new MATRIXEPETRA(M_velocityFESpace->map()));
    SHP(MATRIXEPETRA) jac01(new MATRIXEPETRA(M_velocityFESpace->map()));
    SHP(MATRIXEPETRA) jac10(new MATRIXEPETRA(M_pressureFESpace->map()));
    SHP(MATRIXEPETRA) jac11(new MATRIXEPETRA(M_pressureFESpace->map()));

    jac00->zero();
    jac01->zero();
    jac10->zero();
    jac11->zero();

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
              ) >> jac00;

    integrate(elements(M_velocityFESpace->mesh()),
              M_pressureFESpace->qr(),
              M_pressureFESpaceETA,
              M_velocityFESpaceETA,
              TAU_M * (dot(grad(phi_i),
                       value(M_density) *
                       phi_j * grad(M_velocityFESpaceETA, *velocityRep)
                    +  value(M_density) *
                       value(M_velocityFESpaceETA, *velocityRep) * grad(phi_j)
                    -  value(M_viscosity) *
                       laplacian(phi_j)))
              ) >> jac10;

    integrate(elements(M_velocityFESpace->mesh()),
              M_velocityFESpace->qr(),
              M_velocityFESpaceETA,
              M_pressureFESpaceETA,
              TAU_M * (value(M_density)
                       * dot(value(M_velocityFESpaceETA, *velocityRep)
                       * grad(phi_i), grad(phi_j)))
              ) >> jac01;

    integrate(elements(M_velocityFESpace->mesh()),
              M_pressureFESpace->qr(),
              M_pressureFESpaceETA,
              M_pressureFESpaceETA,
              TAU_M * dot(grad(phi_i), grad(phi_j))
          ) >> jac11;


    jac00->globalAssemble();

    jac10->globalAssemble(M_velocityFESpace->mapPtr(),
                          M_pressureFESpace->mapPtr());

    jac01->globalAssemble(M_pressureFESpace->mapPtr(),
                          M_velocityFESpace->mapPtr());

    jac11->globalAssemble();

    M_jac.resize(2,2);
    M_jac.block(0,0).data() = jac00;
    M_jac.block(0,1).data() = jac01;
    M_jac.block(1,0).data() = jac10;
    M_jac.block(1,1).data() = jac11;

    return M_jac;
}

BlockMatrix<MatrixEp>
SUPGStabilization::
assembleMass(const BlockVector<VectorEp>& sol)
{
    using namespace LifeV;
    using namespace LifeV::ExpressionAssembly;

    SHP(VECTOREPETRA) velocityRep(new VECTOREPETRA(*sol.block(0).data(), Repeated));

    SHP(MATRIXEPETRA) mass;
    mass.reset(new MATRIXEPETRA(M_velocityFESpace->map()));
    mass->zero();
    SHP(SquareRoot) squareroot(new SquareRoot());

    integrate(elements(M_velocityFESpace->mesh()),
              M_velocityFESpace->qr(),
              M_velocityFESpaceETA,
              M_velocityFESpaceETA,
              TAU_M * value(M_density*M_density)
                    * dot(phi_j*grad(phi_i),
                      value(M_velocityFESpaceETA, *velocityRep))
              ) >> mass;
    mass->globalAssemble();

    BlockMatrix<MatrixEp> retMatrix(2,2);
    retMatrix.block(0,0).data() = mass;

    return retMatrix;
}

BlockVector<VectorEp>
SUPGStabilization::
getResidual(const BlockVector<VectorEp>& sol,
            const BlockVector<VectorEp>& rhs)
{
    using namespace LifeV;
    using namespace LifeV::ExpressionAssembly;

    SHP(VECTOREPETRA) velocityRep(new VECTOREPETRA(*sol.block(0).data(), Repeated));
    SHP(VECTOREPETRA) velocityRhsRep(new VECTOREPETRA(*rhs.block(0).data(), Repeated));
    SHP(VECTOREPETRA) pressureRep(new VECTOREPETRA(*sol.block(1).data(), Repeated));

    SHP(SquareRoot) squareroot(new SquareRoot());

    SHP(VECTOREPETRA) resvelrep(new VECTOREPETRA(M_velocityFESpace->map(), Repeated));
    resvelrep->zero();

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
            ) >> resvelrep;
    resvelrep->globalAssemble();
    SHP(VECTOREPETRA) resvel(new VECTOREPETRA(*resvelrep, Unique));

    SHP(VECTOREPETRA) respressrep(new VECTOREPETRA(M_pressureFESpace->map(), Repeated));
    respressrep->zero();

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
             ) >> respressrep;
    respressrep->globalAssemble();
    SHP(VECTOREPETRA) respress(new VECTOREPETRA(*respressrep, Unique));

    BlockVector<VectorEp> retVec(2);
    retVec.block(0).data() = resvel;
    retVec.block(1).data() = respress;

    return retVec;
}


}  // namespace RedMA
