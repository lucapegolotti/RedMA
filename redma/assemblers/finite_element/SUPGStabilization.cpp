// implementation inspired by StabilizationSUPG.cpp in LifeV
#include "SUPGStabilization.hpp"

namespace RedMA
{

SUPGStabilization::
SUPGStabilization(const DataContainer& data,
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
SUPGStabilization::
getMass(SHP(BlockVector) sol, SHP(BlockVector) rhs)
{
    using namespace LifeV;
    using namespace LifeV::ExpressionAssembly;
    SHP(SquareRoot) squareroot(new SquareRoot());

    SHP(VECTOREPETRA) velocityRep(new VECTOREPETRA(*std::static_pointer_cast<VECTOREPETRA>(sol->block(0)->data()), Repeated));
    SHP(VECTOREPETRA) velocityRhsRep(new VECTOREPETRA(*std::static_pointer_cast<VECTOREPETRA>(rhs->block(0)->data()), Repeated));
    SHP(VECTOREPETRA) pressureRep(new VECTOREPETRA(*std::static_pointer_cast<VECTOREPETRA>(sol->block(1)->data()), Repeated));

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

    SHP(SparseMatrix) mass00wrap(new SparseMatrix());
    mass00wrap->setMatrix(mass00);

    SHP(SparseMatrix) mass10wrap(new SparseMatrix());
    mass10wrap->setMatrix(mass10);

    SHP(SparseMatrix) mass01wrap(new SparseMatrix());
    mass01wrap->setMatrix(mass01);

    SHP(SparseMatrix) mass11wrap(new SparseMatrix());
    mass11wrap->setMatrix(mass11);

    M_mass.reset(new BlockMatrix(2,2));
    M_mass->setBlock(0,0,mass00wrap);
    M_mass->setBlock(0,1,mass01wrap);
    M_mass->setBlock(1,0,mass10wrap);
    M_mass->setBlock(1,1,mass11wrap);

    return M_mass;
}

SHP(BlockMatrix)
SUPGStabilization::
getMassJac(SHP(BlockVector) sol, SHP(BlockVector) rhs)
{
    using namespace LifeV;
    using namespace LifeV::ExpressionAssembly;
    SHP(SquareRoot) squareroot(new SquareRoot());

    SHP(VECTOREPETRA) velocityRep(new VECTOREPETRA(*std::static_pointer_cast<VECTOREPETRA>(sol->block(0)->data()), Repeated));
    SHP(VECTOREPETRA) velocityRhsRep(new VECTOREPETRA(*std::static_pointer_cast<VECTOREPETRA>(rhs->block(0)->data()), Repeated));
    SHP(VECTOREPETRA) pressureRep(new VECTOREPETRA(*std::static_pointer_cast<VECTOREPETRA>(sol->block(1)->data()), Repeated));

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
                     * dot(phi_j * grad(phi_i),
                       value(VH)))
             ) >> massjac00;

    massjac00->globalAssemble();

    massjac10->globalAssemble(M_velocityFESpace->mapPtr(),
                              M_pressureFESpace->mapPtr());

    massjac01->globalAssemble(M_pressureFESpace->mapPtr(),
                              M_velocityFESpace->mapPtr());

    massjac11->globalAssemble();

    SHP(SparseMatrix) massJac00wrap(new SparseMatrix());
    massJac00wrap->setMatrix(massjac00);

    SHP(SparseMatrix) massJac10wrap(new SparseMatrix());
    massJac10wrap->setMatrix(massjac10);

    SHP(SparseMatrix) massJac01wrap(new SparseMatrix());
    massJac01wrap->setMatrix(massjac01);

    SHP(SparseMatrix) massJac11wrap(new SparseMatrix());
    massJac11wrap->setMatrix(massjac11);

    M_massJac.reset(new BlockMatrix(2,2));
    M_massJac->setBlock(0,0,massJac00wrap);
    M_massJac->setBlock(0,1,massJac01wrap);
    M_massJac->setBlock(1,0,massJac10wrap);
    M_massJac->setBlock(1,1,massJac11wrap);

    return M_massJac;
}

SHP(BlockMatrix)
SUPGStabilization::
getJac(SHP(BlockVector) sol, SHP(BlockVector) rhs)
{
    using namespace LifeV;
    using namespace LifeV::ExpressionAssembly;
    SHP(SquareRoot) squareroot(new SquareRoot());

    SHP(VECTOREPETRA) velocityRep(new VECTOREPETRA(*std::static_pointer_cast<VECTOREPETRA>(sol->block(0)->data()), Repeated));
    SHP(VECTOREPETRA) velocityRhsRep(new VECTOREPETRA(*std::static_pointer_cast<VECTOREPETRA>(rhs->block(0)->data()), Repeated));
    SHP(VECTOREPETRA) pressureRep(new VECTOREPETRA(*std::static_pointer_cast<VECTOREPETRA>(sol->block(1)->data()), Repeated));

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

    SHP(SparseMatrix) jac00wrap(new SparseMatrix());
    jac00wrap->setMatrix(jac00);

    SHP(SparseMatrix) jac10wrap(new SparseMatrix());
    jac10wrap->setMatrix(jac10);

    SHP(SparseMatrix) jac01wrap(new SparseMatrix());
    jac01wrap->setMatrix(jac01);

    SHP(SparseMatrix) jac11wrap(new SparseMatrix());
    jac11wrap->setMatrix(jac11);

    M_jac.reset(new BlockMatrix(2,2));
    M_jac->resize(2,2);
    M_jac->setBlock(0,0,jac00wrap);
    M_jac->setBlock(0,1,jac01wrap);
    M_jac->setBlock(1,0,jac10wrap);
    M_jac->setBlock(1,1,jac11wrap);

    return M_jac;
}

SHP(BlockVector)
SUPGStabilization::
getResidual(SHP(BlockVector) sol,
            SHP(BlockVector) rhs)
{
    using namespace LifeV;
    using namespace LifeV::ExpressionAssembly;

    SHP(VECTOREPETRA) velocityRep(new VECTOREPETRA(*std::static_pointer_cast<VECTOREPETRA>(sol->block(0)->data()), Repeated));
    SHP(VECTOREPETRA) velocityRhsRep(new VECTOREPETRA(*std::static_pointer_cast<VECTOREPETRA>(rhs->block(0)->data()), Repeated));
    SHP(VECTOREPETRA) pressureRep(new VECTOREPETRA(*std::static_pointer_cast<VECTOREPETRA>(sol->block(1)->data()), Repeated));

    SHP(SquareRoot) squareroot(new SquareRoot());

    SHP(VECTOREPETRA) resvelrep(new VECTOREPETRA(M_velocityFESpace->map(), Repeated));
    resvelrep->zero();

    integrate(elements(M_velocityFESpace->mesh()),
              M_velocityFESpace->qr(),
              M_velocityFESpaceETA,
              TAU_M * (dot(value(M_density) * value(VH) * grad(phi_i),
                       MOMENTUM_R))
             +TAU_C * div(phi_i) * trace(grad(VH))
            ) >> resvelrep;
    resvelrep->globalAssemble();
    SHP(VECTOREPETRA) resvel(new VECTOREPETRA(*resvelrep, Unique));

    SHP(VECTOREPETRA) respressrep(new VECTOREPETRA(M_pressureFESpace->map(), Repeated));
    respressrep->zero();

    integrate(elements(M_velocityFESpace->mesh()),
              M_pressureFESpace->qr(),
              M_pressureFESpaceETA,
              TAU_M * (dot(grad(phi_i),MOMENTUM_R))
             ) >> respressrep;
    respressrep->globalAssemble();
    SHP(VECTOREPETRA) respress(new VECTOREPETRA(*respressrep, Unique));

    SHP(DistributedVector) comp0(new DistributedVector());
    comp0->setVector(resvel);

    SHP(DistributedVector) comp1(new DistributedVector());
    comp1->setVector(respress);

    SHP(BlockVector) retVec(new BlockVector(2));
    retVec->setBlock(0,comp0);
    retVec->setBlock(1,comp1);

    return retVec;
}


}  // namespace RedMA
