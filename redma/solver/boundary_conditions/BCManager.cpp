#include "BCManager.hpp"

namespace RedMA
{

BCManager::
BCManager(const DataContainer& data, SHP(TreeNode) treeNode) :
  M_data(data),
  M_treeNode(treeNode)
{
    M_inflow = data.getInflow();
    M_strongDirichlet = std::strcmp(data("bc_conditions/inletdirichlet", "weak").c_str(),"strong") == 0;

    parseNeumannData();
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
addInletBC(SHP(LifeV::BCHandler) bcs, std::function<double(double)> law) const
{
    if (M_treeNode->isInletNode())
    {
        auto foo = std::bind(poiseulle,
                             std::placeholders::_1,
                             std::placeholders::_2,
                             std::placeholders::_3,
                             std::placeholders::_4,
                             std::placeholders::_5,
                             M_treeNode->M_block->getInlet(),
                             law);

        LifeV::BCFunctionBase inflowFunction(foo);
        bcs->addBC("Inlet", inletFlag, LifeV::Essential, LifeV::Full,
                   inflowFunction, 3);
    }
}

void
BCManager::
applyDirichletBCs(const double& time, BlockVector<VectorEp>& input,
                  SHP(FESPACE) fespace, const unsigned int& index) const
{
    SHP(LifeV::BCHandler) bcs = createBCHandler0Dirichlet();

    addInletBC(bcs, M_inflow);

    bcs->bcUpdate(*fespace->mesh(), fespace->feBd(), fespace->dof());

    auto curVec = input.block(index).data();

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
apply0DirichletMatrix(BlockMatrix<MatrixEp>& input,
                      SHP(FESPACE) fespace,
                      const unsigned int& index,
                      const double& diagCoefficient) const
{
    SHP(LifeV::BCHandler) bcs = createBCHandler0Dirichlet();

    if (M_strongDirichlet)
        addInletBC(bcs, fZero2);

    bcs->bcUpdate(*fespace->mesh(), fespace->feBd(), fespace->dof());

    unsigned int nCols = input.nCols();

    for (unsigned int j = 0; j < nCols; j++)
    {
        auto curMatrix = input.block(index, j).data();
        if (curMatrix)
        {
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
apply0DirichletBCs(BlockVector<VectorEp>& input, SHP(FESPACE) fespace,
                   const unsigned int& index) const
{
    SHP(LifeV::BCHandler) bcs = createBCHandler0Dirichlet();

    if (M_strongDirichlet)
        addInletBC(bcs, fZero2);

    bcs->bcUpdate(*fespace->mesh(), fespace->feBd(), fespace->dof());

    auto curVec = input.block(index).data();
    if (curVec)
        bcManageRhs(*curVec, *fespace->mesh(), fespace->dof(),
                    *bcs, fespace->feBd(), 0.0, 0.0);
}

SHP(LifeV::BCHandler)
BCManager::
createBCHandler0Dirichlet() const
{
    LifeV::BCFunctionBase zeroFunction(fZero);

    SHP(LifeV::BCHandler) bcs;
    bcs.reset(new LifeV::BCHandler);

    bcs->addBC("Wall", wallFlag, LifeV::Essential,
               LifeV::Full, zeroFunction, 3);
    bcs->addBC("InletRing", inletRing, LifeV::EssentialEdges,
               LifeV::Full, zeroFunction, 3);
    bcs->addBC("OutletRing", outletRing, LifeV::EssentialEdges,
               LifeV::Full, zeroFunction, 3);

    return bcs;
}

double
BCManager::
poiseulle(const double& t, const double& x, const double& y,
          const double& z, const unsigned int& i,
          const GeometricFace& face, const std::function<double(double)> inflow)
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

    Vector3D inflowValue = -inflowNorm * normal;
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
