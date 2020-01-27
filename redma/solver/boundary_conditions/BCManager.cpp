#include "BCManager.hpp"

namespace RedMA
{

BCManager::
BCManager(const DataContainer& data, SHP(TreeNode) treeNode) :
  M_data(data),
  M_treeNode(treeNode)
{
    M_inflow = data.getInflow();
}

void
BCManager::
applyDirichletBCs(const double& time, BlockVector<VectorEp>& input,
                  SHP(FESPACE) fespace, const unsigned int& index) const
{
    SHP(LifeV::BCHandler) bcs = createBCHandler0Dirichlet();

    if (M_treeNode->M_ID == 0)
    {
        auto foo = std::bind(poiseulle,
                             std::placeholders::_1,
                             std::placeholders::_2,
                             std::placeholders::_3,
                             std::placeholders::_4,
                             std::placeholders::_5,
                             M_treeNode->M_block->getInlet(),
                             M_inflow);

        LifeV::BCFunctionBase inflowFunction(foo);
        bcs->addBC("Inlet", inletFlag, LifeV::Essential, LifeV::Full,
                   inflowFunction, 3);
    }

    bcs->bcUpdate(*fespace->mesh(), fespace->feBd(), fespace->dof());

    bcManageRhs(*input.block(index).data(), *fespace->mesh(), fespace->dof(),
                *bcs, fespace->feBd(), 1.0, time);
}

void
BCManager::
apply0DirichletMatrix(BlockMatrix<MatrixEp>& input,
                      SHP(FESPACE) fespace,
                      const unsigned int& index,
                      const double& diagCoefficient) const
{
    SHP(LifeV::BCHandler) bcs = createBCHandler0Dirichlet();
    bcs->bcUpdate(*fespace->mesh(), fespace->feBd(), fespace->dof());

    unsigned int nCols = input.nCols();

    for (unsigned int j = 0; j < nCols; j++)
    {
        auto curMatrix = input.block(index, j).data();
        if (curMatrix)
            bcManageMatrix(*curMatrix, *fespace->mesh(),
                           fespace->dof(), *bcs, fespace->feBd(),
                           (j == index) * diagCoefficient, 0.0);
    }

}

void
BCManager::
apply0DirichletBCs(BlockVector<VectorEp>& input, SHP(FESPACE) fespace,
                   const unsigned int& index) const
{
    SHP(LifeV::BCHandler) bcs = createBCHandler0Dirichlet();

    bcs->bcUpdate(*fespace->mesh(), fespace->feBd(), fespace->dof());

    bcManageRhs(*input.block(index).data(), *fespace->mesh(), fespace->dof(),
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

void
BCManager::
setInflow(std::function<double(double)> inflow)
{
    M_inflow = inflow;
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
    double radius = face.M_radius;

    Vector3D curPoint(x,y,z);
    Vector3D diff = curPoint - center;
    double normDiff = radius - diff.norm();

    double inflowNorm = inflow(t) * normDiff * normDiff / (radius * radius);

    Vector3D inflowValue = -inflowNorm * normal;
    return inflowValue[i];
}

double
BCManager::
fZero(const double& t, const double& x, const double& y,
      const double& z, const unsigned int& i)
{
    return 0.0;
}

}
