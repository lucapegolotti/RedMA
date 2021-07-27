#include "IPStabilization.hpp"

namespace RedMA {

IPStabilization::
IPStabilization(const DataContainer &data,
                shp<FESPACE> fespaceVelocity,
                shp<FESPACE> fespacePressure,
                shp<ETFESPACE3> etfespaceVelocity,
                shp<ETFESPACE1> etfespacePressure,
                EPETRACOMM comm):
    NavierStokesStabilization(data,
                              fespaceVelocity, fespacePressure,
                              etfespaceVelocity, etfespacePressure,
                              comm),
    M_pressureFESpaceSide1(new LifeV::CurrentFE(M_pressureFESpace->refFE(),
                           getGeometricMap(*(M_pressureFESpace->mesh())),
                           M_pressureFESpace->qr())),
    M_pressureFESpaceSide2(new LifeV::CurrentFE(M_pressureFESpace->refFE(),
                           getGeometricMap(*(M_pressureFESpace->mesh())),
                           M_pressureFESpace->qr()))
{
    M_gamma_pressure = data("assembler/stabilization/gamma_pressure", 0.2);

    if (M_gamma_pressure < 0)
        throw new Exception("IP stabilization constants must be positive!");
}

void
IPStabilization::
setup()
{
    getMass();
    getMassJacobian();
    getJacobian();
}

shp<BlockMatrix>
IPStabilization::
getMass(shp<BlockVector> sol,
        shp<BlockVector> rhs)
{
    if (!M_mass)
        M_mass.reset(new BlockMatrix(2,2));

    return M_mass;
}

shp<BlockMatrix>
IPStabilization::
getMassJacobian(shp<BlockVector> sol,
                shp<BlockVector> rhs)
{
    if (!M_massJac)
        M_massJac.reset(new BlockMatrix(2,2));

    return M_massJac;
}

shp<BlockMatrix>
IPStabilization::
getJacobian(shp<BlockVector> sol,
            shp<BlockVector> rhs)
{
    using namespace LifeV;

    if (M_jac)
        return M_jac;

    shp<MATRIXEPETRA> jac11;
    jac11.reset(new MATRIXEPETRA(*(this->stabilizePressure())));

    M_jac.reset(new BlockMatrix(2,2));
    M_jac->setBlock(1,1, wrap(jac11));

    return M_jac;
}

shp<BlockVector>
IPStabilization::
getResidual(shp<BlockVector> sol,
            shp<BlockVector> rhs)
{
    if (!M_jac)
        getJacobian();

    shp<BlockVector> retVec = convert<BlockVector>(M_jac->multiplyByVector(sol));

    return retVec;
}

/*shp<MATRIXEPETRA>
IPStabilization::
stabilizePressure()
{
    using namespace LifeV;
    using namespace LifeV::ExpressionAssembly;

    shp<MATRIXEPETRA> A(new MATRIXEPETRA(M_pressureFESpace->map()));

    integrate(elements(M_velocityFESpaceETA->mesh()),
              M_velocityFESpace->qr(),
              M_pressureFESpaceETA,
              M_pressureFESpaceETA,
              dot(grad(phi_i), grad(phi_j))
    ) >> A;

    A->globalAssemble();

    return A;
}*/

shp<MATRIXEPETRA>
IPStabilization::
stabilizePressure()
{
    using namespace LifeV;

    using faceList_Type = std::vector<MESH::facet_Type*>;

    shp<MATRIXEPETRA> jac(new MATRIXEPETRA(M_pressureFESpace->map()));
    jac->zero();

    std::shared_ptr<faceList_Type> faceInteriorListPtr(new faceList_Type());
    for (auto & face : M_pressureFESpace->mesh()->faceList)
        if (!face.boundary())
            faceInteriorListPtr->push_back(&face);

    // ID geoDims = MESH::S_geoDimensions;
    // const unsigned int nDof = M_pressureFESpace->dofPtr()->numTotalDof();

    unsigned int myFacets = 0;
    for (auto & face : *faceInteriorListPtr)
    {
        const unsigned int iElAd1(face->firstAdjacentElementIdentity());
        const unsigned int iElAd2(face->secondAdjacentElementIdentity());

        if (!(myFacets % 1000)) {
            std::string msg = "[IPStabilization] Stabilizing pressure at internal facet " +
                              std::to_string(myFacets+1) + " out of " +
                              std::to_string(faceInteriorListPtr->size())
                              + "...\n";
            printlog(YELLOW, msg, true);
        }
        ++myFacets;

        (M_pressureFESpace->feBd()).update(*face, UPDATE_W_ROOT_DET_METRIC | UPDATE_QUAD_NODES);

        M_pressureFESpaceSide1->updateFirstDeriv(M_pressureFESpace->mesh()->element(iElAd1));
        M_pressureFESpaceSide2->updateFirstDeriv(M_pressureFESpace->mesh()->element(iElAd2));

        MatrixElemental elMatP(M_pressureFESpace->dof().numTotalDof(),
                               1, 1);

        const double hK2 = M_pressureFESpace->feBd().measure();
        double coeff_pressure = M_gamma_pressure * M_density / M_viscosity *
                                hK2 * std::sqrt(hK2);

        // coeff_press * int_{Face} grad(p)_1 . grad(q)_1
        elMatP.zero();
        AssemblyElemental::ipstab_grad(coeff_pressure, elMatP,
                                       *M_pressureFESpaceSide1, *M_pressureFESpaceSide1,
                                       M_pressureFESpace->feBd(),
                                       0, 0);
        assembleMatrix(*jac, elMatP,
                       *M_pressureFESpaceSide1, *(M_pressureFESpace->dofPtr()),
                       0, 0, 0, 0);

        // coeff_press * int_{Face} grad(p)_1 . grad(q)_1
        elMatP.zero();
        AssemblyElemental::ipstab_grad(coeff_pressure, elMatP,
                                       *M_pressureFESpaceSide2, *M_pressureFESpaceSide2,
                                       M_pressureFESpace->feBd(),
                                       0, 0);
        assembleMatrix(*jac, elMatP,
                       *M_pressureFESpaceSide2, *(M_pressureFESpace->dofPtr()),
                       0, 0, 0, 0);

        // - coeff_press * int_{Face} grad(p)_1 . grad(q)_2
        elMatP.zero();
        AssemblyElemental::ipstab_grad(-coeff_pressure, elMatP,
                                       *M_pressureFESpaceSide1, *M_pressureFESpaceSide2,
                                       M_pressureFESpace->feBd(),
                                       0, 0);
        assembleMatrix(*jac, elMatP,
                       *M_pressureFESpaceSide1, *M_pressureFESpaceSide2,
                       *(M_pressureFESpace->dofPtr()), *(M_pressureFESpace->dofPtr()),
                       0, 0, 0, 0);

        // - coeff_press * int_{Face} grad(p)_2 . grad(q)_1
        elMatP.zero();
        AssemblyElemental::ipstab_grad(-coeff_pressure, elMatP,
                                       *M_pressureFESpaceSide2, *M_pressureFESpaceSide1,
                                       M_pressureFESpace->feBd(),
                                       0, 0);
        assembleMatrix(*jac, elMatP,
                       *M_pressureFESpaceSide2, *M_pressureFESpaceSide1,
                       *(M_pressureFESpace->dofPtr()), *(M_pressureFESpace->dofPtr()),
                       0, 0, 0, 0);
    }

    jac->globalAssemble();

    return jac;
}

}  // namespace RedMA
