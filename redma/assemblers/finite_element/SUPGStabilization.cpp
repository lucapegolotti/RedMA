// implementation inspired by StabilizationSUPG.cpp in LifeV
#include "SUPGStabilization.hpp"

namespace RedMA
{

SUPGStabilization::
SUPGStabilization(const DataContainer& data,
                  shp<FESPACE> fespaceVelocity,
                  shp<FESPACE> fespacePressure,
                  shp<ETFESPACE3> etfespaceVelocity,
                  shp<ETFESPACE1> etfespacePressure) :
  NavierStokesStabilization(data,
                            fespaceVelocity, fespacePressure,
                            etfespaceVelocity, etfespacePressure)
{
}

shp<BlockMatrix>
SUPGStabilization::
getMass(shp<BlockVector> sol, shp<BlockVector> rhs)
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

    shp<SparseMatrix> mass00wrap(new SparseMatrix());
    mass00wrap->setMatrix(mass00);

    shp<SparseMatrix> mass10wrap(new SparseMatrix());
    mass10wrap->setMatrix(mass10);

    shp<SparseMatrix> mass01wrap(new SparseMatrix());
    mass01wrap->setMatrix(mass01);

    shp<SparseMatrix> mass11wrap(new SparseMatrix());
    mass11wrap->setMatrix(mass11);

    M_mass.reset(new BlockMatrix(2,2));
    M_mass->setBlock(0,0,mass00wrap);
    M_mass->setBlock(0,1,mass01wrap);
    M_mass->setBlock(1,0,mass10wrap);
    M_mass->setBlock(1,1,mass11wrap);

    return M_mass;
}

shp<BlockMatrix>
SUPGStabilization::
getMassJac(shp<BlockVector> sol, shp<BlockVector> rhs)
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

    shp<SparseMatrix> massJac00wrap(new SparseMatrix());
    massJac00wrap->setMatrix(massjac00);

    shp<SparseMatrix> massJac10wrap(new SparseMatrix());
    massJac10wrap->setMatrix(massjac10);

    shp<SparseMatrix> massJac01wrap(new SparseMatrix());
    massJac01wrap->setMatrix(massjac01);

    shp<SparseMatrix> massJac11wrap(new SparseMatrix());
    massJac11wrap->setMatrix(massjac11);

    M_massJac.reset(new BlockMatrix(2,2));
    M_massJac->setBlock(0,0,massJac00wrap);
    M_massJac->setBlock(0,1,massJac01wrap);
    M_massJac->setBlock(1,0,massJac10wrap);
    M_massJac->setBlock(1,1,massJac11wrap);

    return M_massJac;
}

shp<BlockMatrix>
SUPGStabilization::
getJac(shp<BlockVector> sol, shp<BlockVector> rhs)
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

    shp<SparseMatrix> jac00wrap(new SparseMatrix());
    jac00wrap->setMatrix(jac00);

    shp<SparseMatrix> jac10wrap(new SparseMatrix());
    jac10wrap->setMatrix(jac10);

    shp<SparseMatrix> jac01wrap(new SparseMatrix());
    jac01wrap->setMatrix(jac01);

    shp<SparseMatrix> jac11wrap(new SparseMatrix());
    jac11wrap->setMatrix(jac11);

    M_jac.reset(new BlockMatrix(2,2));
    M_jac->resize(2,2);
    M_jac->setBlock(0,0,jac00wrap);
    M_jac->setBlock(0,1,jac01wrap);
    M_jac->setBlock(1,0,jac10wrap);
    M_jac->setBlock(1,1,jac11wrap);

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

    shp<DistributedVector> comp0(new DistributedVector());
    comp0->setVector(resvel);

    shp<DistributedVector> comp1(new DistributedVector());
    comp1->setVector(respress);

    shp<BlockVector> retVec(new BlockVector(2));
    retVec->setBlock(0,comp0);
    retVec->setBlock(1,comp1);

    return retVec;
}


}  // namespace RedMA
