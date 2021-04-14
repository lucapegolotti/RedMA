#include "BCManager.hpp"

namespace RedMA
{

BCManager::
BCManager(const DataContainer& data, shp<TreeNode> treeNode) :
  M_data(data),
  M_treeNode(treeNode)
{
    M_inflows = M_data.getInflows();
    M_coefficientInflow = data("bc_conditions/coefficientinflow", 1.0);
    parseNeumannData();

    M_wallFlag = treeNode->M_block->wallFlag();
    M_inletRing = 30;
    M_outletRing = 31;
}

void
BCManager::
parseNeumannData()
{
    // unsigned int numConditions = M_data("bc_conditions/numoutletbcs", 0);
    //
    // for (unsigned int outletIndex = 0; outletIndex < numConditions; outletIndex++)
    // {
    //     std::string dataEntry = "bc_conditions/outlet" + std::to_string(outletIndex);
    //
    //     unsigned int blockindex = M_data(dataEntry + "/blockindex", 0);
    //     if (M_treeNode->M_ID == blockindex)
    //     {
    //         unsigned int boundaryflag = M_data(dataEntry + "/boundaryflag", 2);
    //         M_models[boundaryflag].reset(new WindkesselModel(M_data, dataEntry, outletIndex));
    //     }
    // }
}

void
BCManager::
addInletBC(shp<LifeV::BCHandler> bcs, const Law& law,
           GeometricFace inlet) const
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
                             inlet,
                             law,
                             M_coefficientInflow);

        LifeV::BCFunctionBase inflowFunction(foo);
        bcs->addBC("Inlet", inlet.M_flag, LifeV::Essential, LifeV::Full,
                   inflowFunction, 3);
        bcs->addBC("Wall", M_wallFlag, LifeV::Essential,
                   LifeV::Full, zeroFunction, 3);
    }
}

void
BCManager::
applyDirichletBCs(const double& time, BlockVector& input,
                  shp<FESPACE> fespace, const unsigned int& index) const
{
    shp<LifeV::BCHandler> bcs = createBCHandler0Dirichlet();

    if (M_inflows.find(0) != M_inflows.end())
    {
        Law inflow = M_inflows.find(0)->second;
        addInletBC(bcs, inflow, M_treeNode->M_block->getInlet(0));
    }
    else
    {
        auto inlets = M_treeNode->M_block->getInlets();
        for (auto inlet : inlets)
            addInletBC(bcs, M_inflows.find(inlet.M_flag)->second, inlet);
    }

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
                      const double& diagCoefficient) const
{
    shp<LifeV::BCHandler> bcs = createBCHandler0Dirichlet();

    // if (M_strongDirichlet)
    //     addInletBC(bcs, fZero2);

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
                   const unsigned int& index) const
{
    shp<LifeV::BCHandler> bcs = createBCHandler0Dirichlet();

    // if (M_strongDirichlet)
    //     addInletBC(bcs, fZero2);

    bcs->bcUpdate(*fespace->mesh(), fespace->feBd(), fespace->dof());

    shp<VECTOREPETRA> curVec(spcast<VECTOREPETRA>(input.block(index)->data()));
    if (curVec)
        bcManageRhs(*curVec, *fespace->mesh(), fespace->dof(),
                    *bcs, fespace->feBd(), 0.0, 0.0);
}

shp<LifeV::BCHandler>
BCManager::
createBCHandler0Dirichlet() const
{
    LifeV::BCFunctionBase zeroFunction(fZero);

    shp<LifeV::BCHandler> bcs;
    bcs.reset(new LifeV::BCHandler);

    bcs->addBC("Wall", M_wallFlag, LifeV::Essential,
               LifeV::Full, zeroFunction, 3);
    // bcs->addBC("InletRing", M_wallFlag, LifeV::EssentialEdges,
    //            LifeV::Full, zeroFunction, 3);
    // bcs->addBC("OutletRing", M_wallFlag, LifeV::EssentialEdges,
    //            LifeV::Full, zeroFunction, 3);

    return bcs;
}

double
BCManager::
poiseuille(const double& t, const double& x, const double& y,
          const double& z, const unsigned int& i,
          const GeometricFace& face, const Law inflow,
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
    // auto it = M_models.find(flag);
    // if (it == M_models.end())
    //     return 0.0;
    //
    // return -M_models[flag]->getNeumannCondition(time, rate);
}

double
BCManager::
getNeumannJacobian(const double& time, const double& flag, const double& rate)
{
    // auto it = M_models.find(flag);
    // if (it == M_models.end())
    //     return 0.0;
    //
    // return -M_models[flag]->getNeumannJacobian(time, rate);
}

void
BCManager::
postProcess()
{
    // for (auto windkessel : M_models)
    //     windkessel.second->shiftSolutions();
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
constantFunction(const double& t, const double& x, const double& y,
                 const double& z, const unsigned int& i, const double& K)
{
    return K;
}

}
