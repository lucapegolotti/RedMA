// implementation inspired by StabilizationSUPG.cpp in LifeV
#include "VMSStabilization.hpp"

namespace RedMA
{

VMSStabilization::
VMSStabilization(const DataContainer& data,
                 SHP(FESPACE) fespaceVelocity,
                 SHP(FESPACE) fespacePressure,
                 SHP(ETFESPACE3) etfespaceVelocity,
                 SHP(ETFESPACE1) etfespacePressure) :
  NavierStokesStabilization(data,
                            fespaceVelocity, fespacePressure,
                            etfespaceVelocity, etfespacePressure)
{
}

SHP(BlockMatrix)
VMSStabilization::
getMass(SHP(BlockVector) sol, SHP(BlockVector) rhs)
{
    // using namespace LifeV;
    // using namespace LifeV::ExpressionAssembly;
    // SHP(SquareRoot) squareroot(new SquareRoot());
    //
    // SHP(VECTOREPETRA) velocityRep(new VECTOREPETRA(*sol.block(0).data(), Repeated));
    // SHP(VECTOREPETRA) velocityRhsRep(new VECTOREPETRA(*rhs.block(0).data(), Repeated));
    // SHP(VECTOREPETRA) pressureRep(new VECTOREPETRA(*sol.block(1).data(), Repeated));
    //
    // SHP(MATRIXEPETRA) mass00(new MATRIXEPETRA(M_velocityFESpace->map()));
    // SHP(MATRIXEPETRA) mass01(new MATRIXEPETRA(M_velocityFESpace->map()));
    // SHP(MATRIXEPETRA) mass10(new MATRIXEPETRA(M_pressureFESpace->map()));
    // SHP(MATRIXEPETRA) mass11(new MATRIXEPETRA(M_pressureFESpace->map()));
    //
    // mass00->zero();
    // mass01->zero();
    // mass10->zero();
    // mass11->zero();
    //
    // integrate(elements(M_velocityFESpace->mesh()),
    //           M_velocityFESpace->qr(),
    //           M_velocityFESpaceETA,
    //           M_velocityFESpaceETA,
    //           TAU_M * value(M_density*M_density)
    //                 * dot(value(VH)
    //                 * grad(phi_i), phi_j)
    //           -TAU_M * value(M_density*M_density)
    //                 * dot(phi_j * grad(VH), phi_i)
    //           +TAU_M * TAU_M * value(M_density*M_density*M_density)
    //                 * dot(value(VH) * grad(phi_i), phi_j)
    //           ) >> mass00;
    //
    // // note that in Davide's thesis, the stabilization term for the pressure
    // // has no 1/density in front of the test whereas on wikipedia there is
    // // https://en.wikipedia.org/wiki/Streamline_upwind_Petrov%E2%80%93Galerkin_pressure-stabilizing_Petrov%E2%80%93Galerkin_formulation_for_incompressible_Navier%E2%80%93Stokes_equations
    //
    // integrate(elements(M_velocityFESpace->mesh()),
    //           M_pressureFESpace->qr(),
    //           M_pressureFESpaceETA,
    //           M_velocityFESpaceETA,
    //           TAU_M * value(M_density) * dot(grad(phi_i), phi_j)
    //           ) >> mass10;
    //
    // mass00->globalAssemble();
    //
    // mass01->globalAssemble(M_pressureFESpace->mapPtr(),
    //                        M_velocityFESpace->mapPtr());
    //
    // mass10->globalAssemble(M_velocityFESpace->mapPtr(),
    //                        M_pressureFESpace->mapPtr());
    //
    // mass11->globalAssemble();
    //
    // M_mass.resize(2,2);
    // M_mass.block(0,0).data() = mass00;
    // M_mass.block(0,1).data() = mass01;
    // M_mass.block(1,0).data() = mass10;
    // M_mass.block(1,1).data() = mass11;
    //
    // return M_mass;
}

SHP(BlockMatrix)
VMSStabilization::
getMassJac(SHP(BlockVector) sol, SHP(BlockVector) rhs)
{
    // using namespace LifeV;
    // using namespace LifeV::ExpressionAssembly;
    // SHP(SquareRoot) squareroot(new SquareRoot());
    //
    // SHP(VECTOREPETRA) velocityRep(new VECTOREPETRA(*sol.block(0).data(), Repeated));
    // SHP(VECTOREPETRA) velocityRhsRep(new VECTOREPETRA(*rhs.block(0).data(), Repeated));
    // SHP(VECTOREPETRA) pressureRep(new VECTOREPETRA(*sol.block(1).data(), Repeated));
    //
    // SHP(MATRIXEPETRA) massjac00(new MATRIXEPETRA(M_velocityFESpace->map()));
    // SHP(MATRIXEPETRA) massjac01(new MATRIXEPETRA(M_velocityFESpace->map()));
    // SHP(MATRIXEPETRA) massjac10(new MATRIXEPETRA(M_pressureFESpace->map()));
    // SHP(MATRIXEPETRA) massjac11(new MATRIXEPETRA(M_pressureFESpace->map()));
    //
    // massjac00->zero();
    // massjac01->zero();
    // massjac10->zero();
    // massjac11->zero();
    //
    // integrate(elements(M_velocityFESpace->mesh()),
    //           M_velocityFESpace->qr(),
    //           M_velocityFESpaceETA,
    //           M_velocityFESpaceETA,
    //           TAU_M * (value(M_density*M_density)
    //                  * dot(phi_j * grad(phi_i),
    //                    value(VH)))
    //          -TAU_M * value(M_density*M_density)
    //                 * dot(value(VH) * grad(phi_j), phi_i)
    //          +TAU_M * TAU_M * value(M_density*M_density*M_density)
    //              * dot(phi_j * grad(phi_i), value(VH))
    //          ) >> massjac00;
    //
    // massjac00->globalAssemble();
    //
    // massjac10->globalAssemble(M_velocityFESpace->mapPtr(),
    //                           M_pressureFESpace->mapPtr());
    //
    // massjac01->globalAssemble(M_pressureFESpace->mapPtr(),
    //                           M_velocityFESpace->mapPtr());
    //
    // massjac11->globalAssemble();
    //
    // M_massJac.resize(2,2);
    // M_massJac.block(0,0).data() = massjac00;
    // M_massJac.block(0,1).data() = massjac01;
    // M_massJac.block(1,0).data() = massjac10;
    // M_massJac.block(1,1).data() = massjac11;
    //
    // return M_massJac;
}

SHP(BlockMatrix)
VMSStabilization::
getJac(SHP(BlockVector) sol, SHP(BlockVector) rhs)
{
    // using namespace LifeV;
    // using namespace LifeV::ExpressionAssembly;
    // SHP(SquareRoot) squareroot(new SquareRoot());
    //
    // SHP(VECTOREPETRA) velocityRep(new VECTOREPETRA(*sol.block(0).data(), Repeated));
    // SHP(VECTOREPETRA) velocityRhsRep(new VECTOREPETRA(*rhs.block(0).data(), Repeated));
    // SHP(VECTOREPETRA) pressureRep(new VECTOREPETRA(*sol.block(1).data(), Repeated));
    //
    // SHP(MATRIXEPETRA) jac00(new MATRIXEPETRA(M_velocityFESpace->map()));
    // SHP(MATRIXEPETRA) jac01(new MATRIXEPETRA(M_velocityFESpace->map()));
    // SHP(MATRIXEPETRA) jac10(new MATRIXEPETRA(M_pressureFESpace->map()));
    // SHP(MATRIXEPETRA) jac11(new MATRIXEPETRA(M_pressureFESpace->map()));
    //
    // jac00->zero();
    // jac01->zero();
    // jac10->zero();
    // jac11->zero();
    //
    // integrate(elements(M_velocityFESpace->mesh()),
    //           M_velocityFESpace->qr(),
    //           M_velocityFESpaceETA,
    //           M_velocityFESpaceETA,
    //           TAU_M * (dot(value(M_density) * value(VH) * grad(phi_i),
    //                    MOMENTUM_R_DER)
    //                 +  dot(value(M_density) * phi_j * grad(phi_i), MOMENTUM_R))
    //         - TAU_M * (dot(value(M_density) * MOMENTUM_R * grad(phi_j), phi_i)
    //                 +  dot(value(M_density) * MOMENTUM_R_DER * grad(VH), phi_i))
    //         + TAU_M * TAU_M * (dot(value(M_density) * MOMENTUM_R_DER * grad(phi_i), MOMENTUM_R) +
    //                            dot(value(M_density) * MOMENTUM_R * grad(phi_i), MOMENTUM_R_DER))
    //         + TAU_C *  div(phi_i)*div(phi_j)
    //           ) >> jac00;
    //
    // integrate(elements(M_velocityFESpace->mesh()),
    //           M_pressureFESpace->qr(),
    //           M_pressureFESpaceETA,
    //           M_velocityFESpaceETA,
    //           TAU_M * (dot(grad(phi_i), MOMENTUM_R_DER))
    //           ) >> jac10;
    //
    // integrate(elements(M_velocityFESpace->mesh()),
    //           M_velocityFESpace->qr(),
    //           M_velocityFESpaceETA,
    //           M_pressureFESpaceETA,
    //           TAU_M * (value(M_density) * dot(value(VH) * grad(phi_i), grad(phi_j)))
    //          -TAU_M * (dot(value(M_density) * grad(phi_j) * grad(VH), phi_i))
    //          // there should probably a term like this in the jacobian but it seems to
    //          // lead to worse convergence. Probably it is miscomputed
    //          // +TAU_M * TAU_M * (dot(value(M_density) * MOMENTUM_R * grad(phi_i), grad(phi_j)) +
    //          //                   dot(value(M_density) * grad(phi_j) * grad(phi_i), MOMENTUM_R)
    //          //                   )
    //           ) >> jac01;
    //
    // integrate(elements(M_velocityFESpace->mesh()),
    //           M_pressureFESpace->qr(),
    //           M_pressureFESpaceETA,
    //           M_pressureFESpaceETA,
    //           TAU_M * dot(grad(phi_i), grad(phi_j))
    //       ) >> jac11;
    //
    //
    // jac00->globalAssemble();
    //
    // jac10->globalAssemble(M_velocityFESpace->mapPtr(),
    //                       M_pressureFESpace->mapPtr());
    //
    // jac01->globalAssemble(M_pressureFESpace->mapPtr(),
    //                       M_velocityFESpace->mapPtr());
    //
    // jac11->globalAssemble();
    //
    // M_jac.resize(2,2);
    // M_jac.block(0,0).data() = jac00;
    // M_jac.block(0,1).data() = jac01;
    // M_jac.block(1,0).data() = jac10;
    // M_jac.block(1,1).data() = jac11;
    //
    // return M_jac;
}

SHP(BlockVector)
VMSStabilization::
getResidual(SHP(BlockVector) sol,
            SHP(BlockVector) rhs)
{
    // using namespace LifeV;
    // using namespace LifeV::ExpressionAssembly;
    //
    // SHP(VECTOREPETRA) velocityRep(new VECTOREPETRA(*sol.block(0).data(), Repeated));
    // SHP(VECTOREPETRA) velocityRhsRep(new VECTOREPETRA(*rhs.block(0).data(), Repeated));
    // SHP(VECTOREPETRA) pressureRep(new VECTOREPETRA(*sol.block(1).data(), Repeated));
    //
    // SHP(SquareRoot) squareroot(new SquareRoot());
    //
    // SHP(VECTOREPETRA) resvelrep(new VECTOREPETRA(M_velocityFESpace->map(), Repeated));
    // resvelrep->zero();
    //
    // integrate(elements(M_velocityFESpace->mesh()),
    //           M_velocityFESpace->qr(),
    //           M_velocityFESpaceETA,
    //           TAU_M * (dot(value(M_density) * value(VH) * grad(phi_i),
    //                    MOMENTUM_R))
    //          -TAU_M * (dot(value(M_density) * MOMENTUM_R * grad(VH),
    //                    phi_i))
    //          +TAU_M * TAU_M * (dot(value(M_density) * MOMENTUM_R * grad(phi_i),
    //                            MOMENTUM_R))
    //          +TAU_C * div(phi_i) * trace(grad(VH))
    //         ) >> resvelrep;
    // resvelrep->globalAssemble();
    // SHP(VECTOREPETRA) resvel(new VECTOREPETRA(*resvelrep, Unique));
    //
    // SHP(VECTOREPETRA) respressrep(new VECTOREPETRA(M_pressureFESpace->map(), Repeated));
    // respressrep->zero();
    //
    // integrate(elements(M_velocityFESpace->mesh()),
    //           M_pressureFESpace->qr(),
    //           M_pressureFESpaceETA,
    //           TAU_M * (dot(grad(phi_i),MOMENTUM_R))
    //          ) >> respressrep;
    // respressrep->globalAssemble();
    // SHP(VECTOREPETRA) respress(new VECTOREPETRA(*respressrep, Unique));
    //
    // BlockVector retVec(2);
    // retVec.block(0).data() = resvel;
    // retVec.block(1).data() = respress;
    //
    // return retVec;
}

}  // namespace RedMA
