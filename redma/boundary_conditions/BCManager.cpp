#include "BCManager.hpp"

namespace RedMA
{

BCManager::
BCManager(const DataContainer& data, shp<TreeNode> treeNode) :
  M_data(data),
  M_treeNode(treeNode)
{
    M_inflow = data.getInflow();
    M_strongDirichlet = std::strcmp(data("bc_conditions/inletdirichlet", "weak").c_str(),"strong") == 0;
    M_coefficientInflow = data("bc_conditions/coefficientinflow", 1.0);
    parseNeumannData();
    M_inletFlag = treeNode->M_block->getInlet().M_flag;
    M_wallFlag = treeNode->M_block->wallFlag();
    M_inletRingFlag = treeNode->M_block->getInlet().M_ringFlag;
    M_outletRingFlag = treeNode->M_block->getOutlet(0).M_ringFlag;  // assuming all outlets are the same
}

void
BCManager::
parseNeumannData()
{
    unsigned int numConditions = M_data("bc_conditions/numoutletbcs", 0);

    for (unsigned int outletIndex = 0; outletIndex < numConditions; outletIndex++)
    {
        std::string dataEntry = "bc_conditions/outlet" + std::to_string(outletIndex);

        unsigned int blockindex = M_data(dataEntry + "/blockindex", 0);
        if (M_treeNode->M_ID == blockindex)
        {
            unsigned int boundaryflag = M_data(dataEntry + "/boundaryflag", 2);
            M_models[boundaryflag].reset(new WindkesselModel(M_data, dataEntry, outletIndex));
        }
    }
}

void
BCManager::
addInletBC(shp<LifeV::BCHandler> bcs, std::function<double(double)> law, const bool& ringOnly) const
{
    if (M_treeNode->isInletNode())
    {
        LifeV::BCFunctionBase zeroFunction(fZero);

        auto foo = std::bind(poiseuille,
                             std::placeholders::_1,
                             std::placeholders::_2,
                             std::placeholders::_3,
                             std::placeholders::_4,
                             std::placeholders::_5,
                             M_treeNode->M_block->getInlet(),
                             law,
                             M_coefficientInflow);

        LifeV::BCFunctionBase inflowFunction(foo);
        bcs->addBC("Inlet", M_inletFlag, LifeV::Essential, LifeV::Full,
                   inflowFunction, 3);

        if (ringOnly) {
            bcs->addBC("InletRing", M_inletRingFlag, LifeV::EssentialEdges,
                       LifeV::Full, zeroFunction, 3);
        }
        else {
            bcs->addBC("Wall", M_wallFlag, LifeV::Essential,
                       LifeV::Full, zeroFunction, 3);
        }
    }
}

void
BCManager::
applyDirichletBCs(const double& time, BlockVector& input,
                  shp<FESPACE> fespace, const unsigned int& index,
                  const bool& ringOnly) const
{
    shp<LifeV::BCHandler> bcs = createBCHandler0Dirichlet(ringOnly);

    addInletBC(bcs, M_inflow, ringOnly);

    bcs->bcUpdate(*fespace->mesh(), fespace->feBd(), fespace->dof());

    shp<VECTOREPETRA> curVec(spcast<VECTOREPETRA>(input.block(index)->data()));

    if (curVec)
        bcManageRhs(*curVec, *fespace->mesh(), fespace->dof(),
                    *bcs, fespace->feBd(), 1.0, time);
}

double
BCManager::
fZero2(double time)
{
    return 0.0;
}

void
BCManager::
apply0DirichletMatrix(BlockMatrix& input,
                      shp<FESPACE> fespace,
                      const unsigned int& index,
                      const double& diagCoefficient,
                      const bool& ringOnly) const
{
    shp<LifeV::BCHandler> bcs = createBCHandler0Dirichlet(ringOnly);

    if (M_strongDirichlet)
        addInletBC(bcs, fZero2, ringOnly);

    bcs->bcUpdate(*fespace->mesh(), fespace->feBd(), fespace->dof());

    unsigned int nCols = input.nCols();

    for (unsigned int j = 0; j < nCols; j++)
    {
        if (!input.block(index, j)->isZero())
        {
            shp<MATRIXEPETRA> curMatrix = spcast<MATRIXEPETRA>(input.block(index, j)->data());
            auto domainMap = curMatrix->domainMapPtr();
            auto rangeMap = curMatrix->rangeMapPtr();
            bcManageMatrix(*curMatrix, *fespace->mesh(),
                           fespace->dof(), *bcs, fespace->feBd(),
                           (j == index) * diagCoefficient, 0.0);
            curMatrix->globalAssemble(domainMap, rangeMap);
        }
    }
}

void
BCManager::
apply0DirichletBCs(BlockVector& input, shp<FESPACE> fespace,
                   const unsigned int& index,
                   const bool& ringOnly) const
{
    shp<LifeV::BCHandler> bcs = createBCHandler0Dirichlet(ringOnly);

    if (M_strongDirichlet)
        addInletBC(bcs, fZero2, ringOnly);

    bcs->bcUpdate(*fespace->mesh(), fespace->feBd(), fespace->dof());

    shp<VECTOREPETRA> curVec(spcast<VECTOREPETRA>(input.block(index)->data()));
    if (curVec)
        bcManageRhs(*curVec, *fespace->mesh(), fespace->dof(),
                    *bcs, fespace->feBd(), 0.0, 0.0);
}

shp<LifeV::BCHandler>
BCManager::
createBCHandler0Dirichlet(const bool& ringOnly) const
{
    LifeV::BCFunctionBase zeroFunction(fZero);

    shp<LifeV::BCHandler> bcs;
    bcs.reset(new LifeV::BCHandler);

    if (!(ringOnly))
        bcs->addBC("Wall", M_wallFlag, LifeV::Essential,
                   LifeV::Full, zeroFunction, 3);
    else {
        if (M_treeNode->isInletNode())
            bcs->addBC("InletRing", M_inletRingFlag, LifeV::EssentialEdges,
                       LifeV::Full, zeroFunction, 3);

        if (M_treeNode->isOutletNode())
            // this imposes the BC at all outlets, as the outlet ring flag is the same
            // for all outlets; thus outlet building blocks MUST be tubes!!
            bcs->addBC("OutletRing", M_outletRingFlag, LifeV::EssentialEdges,
                       LifeV::Full, zeroFunction, 3);
    }

    return bcs;
}

double
BCManager::
poiseuille(const double& t, const double& x, const double& y,
          const double& z, const unsigned int& i,
          const GeometricFace& face, const std::function<double(double)> inflow,
          const double& coefficient)
{
    typedef LifeV::VectorSmall<3>   Vector3D;
    // GeometricFace face = M_treeNode->M_block->getInlet();

    const Vector3D& center = face.M_center;
    const Vector3D& normal = face.M_normal;
    double R = face.M_radius;

    Vector3D curPoint(x,y,z);
    Vector3D diff = curPoint - center;
    double r = diff.norm();

    // we suppose that inflow is the flowrate and we want to find the max velocity
    const double maxU = inflow(t) * 2.0 / (M_PI * R * R);
    double inflowNorm = maxU * (1.0 - (r * r)/(R * R));

    if (inflow(t) < 0)
        inflowNorm = inflowNorm < 0 ? inflowNorm : 0;

    inflowNorm = inflowNorm > 0 ? inflowNorm : 0;

    Vector3D inflowValue = -inflowNorm * normal * coefficient;
    return inflowValue[i];
}

double
BCManager::
getNeumannBc(const double& time, const double& flag, const double& rate)
{
    auto it = M_models.find(flag);
    if (it == M_models.end())
        return 0.0;

    return -M_models[flag]->getNeumannCondition(time, rate);
}

double
BCManager::
getNeumannJacobian(const double& time, const double& flag, const double& rate)
{
    auto it = M_models.find(flag);
    if (it == M_models.end())
        return 0.0;

    return -M_models[flag]->getNeumannJacobian(time, rate);
}

void
BCManager::
postProcess()
{
    for (auto windkessel : M_models)
        windkessel.second->shiftSolutions();
}

shp<VECTOREPETRA>
BCManager::
computeBoundaryIndicator(shp<FESPACE> fespace)
{
    shp<VECTOREPETRA> boundaryIndicator(new VECTOREPETRA(fespace->map()));
    boundaryIndicator->zero();

    LifeV::BCFunctionBase oneFunction(fOne);

    shp<LifeV::BCHandler> bcs;
    bcs.reset(new LifeV::BCHandler);

    bcs->addBC("Wall", M_wallFlag, LifeV::Essential,
               LifeV::Full, oneFunction, 3);

    bcs->bcUpdate(*fespace->mesh(), fespace->feBd(), fespace->dof());

    bcManageRhs(*boundaryIndicator, *fespace->mesh(), fespace->dof(),
                *bcs, fespace->feBd(), 1.0, 0.0);

    return boundaryIndicator;
}

double
BCManager::
fZero(const double& t, const double& x, const double& y,
      const double& z, const unsigned int& i)
{
    return 0.0;
}

double
BCManager::
fOne(const double& t, const double& x, const double& y,
      const double& z, const unsigned int& i)
{
    return 1.0;
}

double
BCManager::
constantFunction(const double& t, const double& x, const double& y,
                 const double& z, const unsigned int& i, const double& K)
{
    return K;
}

}
