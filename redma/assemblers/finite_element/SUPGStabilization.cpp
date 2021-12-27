// implementation inspired by StabilizationSUPG.cpp in LifeV
#include "SUPGStabilization.hpp"

namespace RedMA
{

SUPGStabilization::
SUPGStabilization(const DataContainer& data,
                  shp<FESPACE> fespaceVelocity,
                  shp<FESPACE> fespacePressure,
                  shp<ETFESPACE3> etfespaceVelocity,
                  shp<ETFESPACE1> etfespacePressure,
                  EPETRACOMM comm) :
  NavierStokesStabilization(data,
                            fespaceVelocity, fespacePressure,
                            etfespaceVelocity, etfespacePressure,
                            comm)
{
    M_timeOrder = data("time_discretization/order", 2);
    M_dt = data("time_discretization/dt", 0.01);

    if (!std::strcmp(M_velocityOrder.c_str(),"P1"))
        M_C_I = 30;
    else if (!std::strcmp(M_velocityOrder.c_str(),"P2"))
        M_C_I = 60;
    else if (!std::strcmp(M_velocityOrder.c_str(),"P3"))
        M_C_I = 120;
    else if (!std::strcmp(M_velocityOrder.c_str(),"P4"))
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
setup()
{
    // in SUPG stabilization there is nothing to setup in advance
}

shp<BlockMatrix>
SUPGStabilization::
getMass(shp<BlockVector> sol,
        shp<BlockVector> rhs)
{
    using namespace LifeV;
    using namespace LifeV::ExpressionAssembly;
    shp<SquareRoot> squareroot(new SquareRoot());

    shp<VECTOREPETRA> velocityRep(new VECTOREPETRA(*spcast<VECTOREPETRA>(sol->block(0)->data()), Repeated));
    shp<VECTOREPETRA> velocityRhsRep(new VECTOREPETRA(*spcast<VECTOREPETRA>(rhs->block(0)->data()), Repeated));
    shp<VECTOREPETRA> pressureRep(new VECTOREPETRA(*spcast<VECTOREPETRA>(sol->block(1)->data()), Repeated));

    shp<MATRIXEPETRA> mass00(new MATRIXEPETRA(M_velocityFESpace->map()));
    shp<MATRIXEPETRA> mass01(new MATRIXEPETRA(M_velocityFESpace->map()));
    shp<MATRIXEPETRA> mass10(new MATRIXEPETRA(M_pressureFESpace->map()));
    shp<MATRIXEPETRA> mass11(new MATRIXEPETRA(M_pressureFESpace->map()));

    mass00->zero();
    mass01->zero();
    mass10->zero();
    mass11->zero();

    integrate(elements(M_velocityFESpace->mesh()),
              M_velocityFESpace->qr(),
              M_velocityFESpaceETA,
              M_velocityFESpaceETA,
              TAU_M * value(M_density*M_density)
                    * dot(value(VH)
                    * grad(phi_i), phi_j)
              ) >> mass00;

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

    M_mass.reset(new BlockMatrix(2,2));
    M_mass->setBlock(0,0, wrap(mass00));
    M_mass->setBlock(0,1, wrap(mass01));
    M_mass->setBlock(1,0, wrap(mass10));
    M_mass->setBlock(1,1, wrap(mass11));

    return M_mass;
}

shp<BlockMatrix>
SUPGStabilization::
getMassJacobian(shp<BlockVector> sol,
                shp<BlockVector> rhs)
{
    using namespace LifeV;
    using namespace LifeV::ExpressionAssembly;
    shp<SquareRoot> squareroot(new SquareRoot());

    shp<VECTOREPETRA> velocityRep(new VECTOREPETRA(*spcast<VECTOREPETRA>(sol->block(0)->data()), Repeated));
    shp<VECTOREPETRA> velocityRhsRep(new VECTOREPETRA(*spcast<VECTOREPETRA>(rhs->block(0)->data()), Repeated));
    shp<VECTOREPETRA> pressureRep(new VECTOREPETRA(*spcast<VECTOREPETRA>(sol->block(1)->data()), Repeated));

    shp<MATRIXEPETRA> massjac00(new MATRIXEPETRA(M_velocityFESpace->map()));
    shp<MATRIXEPETRA> massjac01(new MATRIXEPETRA(M_velocityFESpace->map()));
    shp<MATRIXEPETRA> massjac10(new MATRIXEPETRA(M_pressureFESpace->map()));
    shp<MATRIXEPETRA> massjac11(new MATRIXEPETRA(M_pressureFESpace->map()));

    massjac00->zero();
    massjac01->zero();
    massjac10->zero();
    massjac11->zero();

    integrate(elements(M_velocityFESpace->mesh()),
              M_velocityFESpace->qr(),
              M_velocityFESpaceETA,
              M_velocityFESpaceETA,
              TAU_M * (value(M_density*M_density)
                     * dot(phi_j * grad(phi_i),
                       value(VH)))
             ) >> massjac00;

    massjac00->globalAssemble();

    massjac10->globalAssemble(M_velocityFESpace->mapPtr(),
                              M_pressureFESpace->mapPtr());

    massjac01->globalAssemble(M_pressureFESpace->mapPtr(),
                              M_velocityFESpace->mapPtr());

    massjac11->globalAssemble();

    M_massJac.reset(new BlockMatrix(2,2));
    M_massJac->setBlock(0,0, wrap(massjac00));
    M_massJac->setBlock(0,1, wrap(massjac01));
    M_massJac->setBlock(1,0, wrap(massjac10));
    M_massJac->setBlock(1,1, wrap(massjac11));

    return M_massJac;
}

shp<BlockMatrix>
SUPGStabilization::
getJacobian(shp<BlockVector> sol,
            shp<BlockVector> rhs)
{
    using namespace LifeV;
    using namespace LifeV::ExpressionAssembly;
    shp<SquareRoot> squareroot(new SquareRoot());

    shp<VECTOREPETRA> velocityRep(new VECTOREPETRA(*spcast<VECTOREPETRA>(sol->block(0)->data()), Repeated));
    shp<VECTOREPETRA> velocityRhsRep(new VECTOREPETRA(*spcast<VECTOREPETRA>(rhs->block(0)->data()), Repeated));
    shp<VECTOREPETRA> pressureRep(new VECTOREPETRA(*spcast<VECTOREPETRA>(sol->block(1)->data()), Repeated));

    shp<MATRIXEPETRA> jac00(new MATRIXEPETRA(M_velocityFESpace->map()));
    shp<MATRIXEPETRA> jac01(new MATRIXEPETRA(M_velocityFESpace->map()));
    shp<MATRIXEPETRA> jac10(new MATRIXEPETRA(M_pressureFESpace->map()));
    shp<MATRIXEPETRA> jac11(new MATRIXEPETRA(M_pressureFESpace->map()));

    jac00->zero();
    jac01->zero();
    jac10->zero();
    jac11->zero();

    integrate(elements(M_velocityFESpace->mesh()),
              M_velocityFESpace->qr(),
              M_velocityFESpaceETA,
              M_velocityFESpaceETA,
              TAU_M * (dot(value(M_density) * value(VH) * grad(phi_i),
                       MOMENTUM_R_DER)
                    +  dot(value(M_density) * phi_j * grad(phi_i), MOMENTUM_R))
            + TAU_C *  div(phi_i)*div(phi_j)
              ) >> jac00;

    integrate(elements(M_velocityFESpace->mesh()),
              M_pressureFESpace->qr(),
              M_pressureFESpaceETA,
              M_velocityFESpaceETA,
              TAU_M * (dot(grad(phi_i), MOMENTUM_R_DER))
              ) >> jac10;

    integrate(elements(M_velocityFESpace->mesh()),
              M_velocityFESpace->qr(),
              M_velocityFESpaceETA,
              M_pressureFESpaceETA,
              TAU_M * dot(value(M_density) * value(VH) * grad(phi_i), grad(phi_j))
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

    M_jac.reset(new BlockMatrix(2,2));
    M_jac->setBlock(0,0, wrap(jac00));
    M_jac->setBlock(0,1, wrap(jac01));
    M_jac->setBlock(1,0, wrap(jac10));
    M_jac->setBlock(1,1, wrap(jac11));

    return M_jac;
}

shp<BlockVector>
SUPGStabilization::
getResidual(shp<BlockVector> sol,
            shp<BlockVector> rhs)
{
    using namespace LifeV;
    using namespace LifeV::ExpressionAssembly;

    shp<VECTOREPETRA> velocityRep(new VECTOREPETRA(*spcast<VECTOREPETRA>(sol->block(0)->data()), Repeated));
    shp<VECTOREPETRA> velocityRhsRep(new VECTOREPETRA(*spcast<VECTOREPETRA>(rhs->block(0)->data()), Repeated));
    shp<VECTOREPETRA> pressureRep(new VECTOREPETRA(*spcast<VECTOREPETRA>(sol->block(1)->data()), Repeated));

    shp<SquareRoot> squareroot(new SquareRoot());

    shp<VECTOREPETRA> resvelrep(new VECTOREPETRA(M_velocityFESpace->map(), Repeated));
    resvelrep->zero();

    integrate(elements(M_velocityFESpace->mesh()),
              M_velocityFESpace->qr(),
              M_velocityFESpaceETA,
              TAU_M * (dot(value(M_density) * value(VH) * grad(phi_i),
                       MOMENTUM_R))
             +TAU_C * div(phi_i) * trace(grad(VH))
            ) >> resvelrep;
    resvelrep->globalAssemble();
    shp<VECTOREPETRA> resvel(new VECTOREPETRA(*resvelrep, Unique));

    shp<VECTOREPETRA> respressrep(new VECTOREPETRA(M_pressureFESpace->map(), Repeated));
    respressrep->zero();

    integrate(elements(M_velocityFESpace->mesh()),
              M_pressureFESpace->qr(),
              M_pressureFESpaceETA,
              TAU_M * (dot(grad(phi_i),MOMENTUM_R))
             ) >> respressrep;
    respressrep->globalAssemble();
    shp<VECTOREPETRA> respress(new VECTOREPETRA(*respressrep, Unique));

    shp<BlockVector> retVec(new BlockVector(2));
    retVec->setBlock(0, wrap(resvel));
    retVec->setBlock(1, wrap(respress));

    return retVec;
}

}  // namespace RedMA
