#include "InterfaceAssembler.hpp"

namespace RedMA
{

template <>
void
InterfaceAssembler<VectorEp, MatrixEp>::
buildMapLagrange(shp<BasisFunctionFunctor> bfs)
{
    unsigned int nBfs = bfs->getNumBasisFunctions();
    auto aFather = M_interface.M_assemblerFather;
    EPETRACOMM comm = aFather->getComm();
    // this must be changed e.g. for scalar coupling ...
    unsigned int myel = (3 * nBfs) / comm->NumProc();

    // the first process takes care of the remainder
    if (comm->MyPID() == 0)
    {
        myel += (3 * nBfs) % comm->NumProc();
    }
    unsigned int mapSize = 3 * nBfs;
    M_mapLagrange.reset(new LifeV::MapEpetra(mapSize, myel, 0, comm));
}

template <>
std::vector<VectorEp>
InterfaceAssembler<VectorEp, MatrixEp>::
buildStabilizationVectorsVelocity(shp<BasisFunctionFunctor> bfs,
                                  const GeometricFace& face,
                                  SHP(aAssembler<VectorEp COMMA MatrixEp>) assembler) const
{
    using namespace LifeV;
    using namespace ExpressionAssembly;

    QuadratureBoundary boundaryQuadRule(buildTetraBDQR
                                       (*generateQuadratureRule("STRANG10")));

    unsigned int nBasisFunctions = bfs->getNumBasisFunctions();

    std::vector<VectorEp> couplingVectors(3 * nBasisFunctions);

    shp<ETFESPACE3> etfespace = assembler->getETFESpaceCoupling();
    MAPEPETRA map = etfespace->map();
    shp<MESH> mesh = etfespace->mesh();
    unsigned int faceFlag = face.M_flag;

    unsigned int count = 0;
    bool useFullStrain = M_data("fluid/use_strain", true);
    double viscosity = this->M_data("fluid/viscosity", 0.035);

    // this must be changed for scalar problems (e.g. laplacian)
    for (unsigned int dim = 0; dim < 3; dim++)
    {
        LifeV::VectorSmall<3> versor;
        versor[0] = 0.;
        versor[1] = 0.;
        versor[2] = 0.;

        versor[dim] = 1.;

        for (unsigned int i = 0; i < nBasisFunctions; i++)
        {
            shp<VECTOREPETRA> currentMode(new VECTOREPETRA(map, LifeV::Repeated));

            bfs->setIndex(i);
            if (useFullStrain)
            {
                integrate(boundary(mesh, faceFlag),
                          boundaryQuadRule,
                          etfespace,
                          eval(bfs, X) * value(0.5 * viscosity) *
                          dot((grad(phi_i) + transpose(grad(phi_i)))*face.M_normal,versor)
                      ) >> currentMode;
            }
            else
            {
                integrate(boundary(mesh, faceFlag),
                          boundaryQuadRule,
                          etfespace,
                          eval(bfs, X) * value(viscosity) *
                          dot((grad(phi_i))*face.M_normal,versor)
                      ) >> currentMode;
            }
            couplingVectors[count].data() = currentMode;
            count++;
        }
    }

    return couplingVectors;
}

template <>
std::vector<VectorEp>
InterfaceAssembler<VectorEp, MatrixEp>::
buildStabilizationVectorsPressure(shp<BasisFunctionFunctor> bfs,
                                  const GeometricFace& face,
                                  SHP(aAssembler<VectorEp COMMA MatrixEp>) assembler) const
{
    using namespace LifeV;
    using namespace ExpressionAssembly;

    QuadratureBoundary boundaryQuadRule(buildTetraBDQR
                                       (*generateQuadratureRule("STRANG10")));

    unsigned int nBasisFunctions = bfs->getNumBasisFunctions();

    std::vector<VectorEp> couplingVectors(3 * nBasisFunctions);

    shp<ETFESPACE1> etfespace = assembler->getETFESpaceSecondary();
    MAPEPETRA map = etfespace->map();
    shp<MESH> mesh = etfespace->mesh();
    unsigned int faceFlag = face.M_flag;

    unsigned int count = 0;

    // this must be changed for scalar problems (e.g. laplacian)
    for (unsigned int dim = 0; dim < 3; dim++)
    {
        LifeV::VectorSmall<3> versor;
        versor[0] = 0.;
        versor[1] = 0.;
        versor[2] = 0.;

        versor[dim] = 1.;

        for (unsigned int i = 0; i < nBasisFunctions; i++)
        {
            shp<VECTOREPETRA> currentMode(new VECTOREPETRA(map, LifeV::Repeated));

            bfs->setIndex(i);

            integrate(boundary(mesh, faceFlag),
                      boundaryQuadRule,
                      etfespace,
                      eval(bfs, X) * value(-1.0) *
                      dot(phi_i*face.M_normal,versor)
                  ) >> currentMode;

            couplingVectors[count].data() = currentMode;
            count++;
        }
    }

    return couplingVectors;
}

template <>
std::vector<VectorEp>
InterfaceAssembler<VectorEp, MatrixEp>::
buildStabilizationVectorsLagrange() const
{
    SHP(const MAPEPETRA) lagrangeMap;
    if (M_fatherBT.nRows() > 0)
        lagrangeMap = M_fatherBT.block(0,0).data()->domainMapPtr();
    else
        lagrangeMap = M_childBT.block(0,0).data()->domainMapPtr();
    unsigned int ncols = lagrangeMap->mapSize();

    std::vector<VectorEp> retVectors(ncols);

    for (unsigned int i = 0; i < ncols; i++)
    {
        shp<VECTOREPETRA> newvec(new VECTOREPETRA(*lagrangeMap, LifeV::Unique));
        newvec->zero();

        if (lagrangeMap->isOwned(i))
            (*newvec)[i] = 1.0;

        retVectors[i].data() = newvec;
    }
    return retVectors;
}

template <>
void
InterfaceAssembler<VectorEp, MatrixEp>::
buildStabilizationMatrix(shp<AssemblerType> assembler,
                         const GeometricFace& face,
                         BlockMatrix<MatrixEp>& matrix)
{
    BlockMatrix<MatrixEp> stab;
    stab.resize(1, assembler->getNumComponents());

    shp<BasisFunctionFunctor> bfs;
    bfs =  BasisFunctionFactory(M_data.getDatafile(), face);

    std::vector<VectorEp> stabVectorsVelocity;
    stabVectorsVelocity = buildStabilizationVectorsVelocity(bfs, face, assembler);

    std::vector<VectorEp> stabVectorsPressure;
    stabVectorsPressure = buildStabilizationVectorsPressure(bfs, face, assembler);

    stab.block(0,0).shallowCopy(MatrixEp(stabVectorsVelocity).transpose());
    stab.block(0,1).shallowCopy(MatrixEp(stabVectorsPressure).transpose());
    matrix.shallowCopy(stab);

    std::vector<VectorEp> stabVectorsLagrange = buildStabilizationVectorsLagrange();
    M_identity.resize(1,1);

    M_identity.block(0,0).shallowCopy(MatrixEp(stabVectorsLagrange));
}

template <>
void
InterfaceAssembler<VectorEp, MatrixEp>::
buildCouplingMatrices()
{
    unsigned int indexOutlet = M_interface.M_indexOutlet;

    auto asFather = M_interface.M_assemblerFather;

    if (asFather)
    {
        GeometricFace outlet = asFather->getTreeNode()->M_block->getOutlet(indexOutlet);

        buildCouplingMatrices(asFather, outlet, M_fatherBT, M_fatherB);

        if (M_stabilizationCoupling > THRESHOLDSTAB)
        {
            buildStabilizationMatrix(asFather, outlet, M_stabFather);
            // M_stabFather *= 0.5;
        }


    }

    auto asChild = M_interface.M_assemblerChild;
    if (asChild)
    {
        GeometricFace inlet = asChild->getTreeNode()->M_block->getInlet();
        // I invert the normal of the face such that it is the same as the outlet
        inlet.M_normal *= (-1.);

        buildCouplingMatrices(asChild, inlet, M_childBT, M_childB);

        if (M_stabilizationCoupling > THRESHOLDSTAB)
        {
            buildStabilizationMatrix(asChild, inlet, M_stabChild);
            // no need to multiply by -1 as the inlet is already reversed
            // M_stabChild *= (-1.);
        }
        M_childB *= (-1.);
        M_childBT *= (-1.);
    }
}

template <>
BlockVector<VectorEp>
InterfaceAssembler<VectorEp, MatrixEp>::
getZeroVector() const
{
    BlockVector<VectorEp> retVector;
    retVector.resize(1);

    shp<MAPEPETRA> lagrangeMap;
    if (M_fatherBT.nRows() > 0)
        lagrangeMap.reset(new MAPEPETRA(*M_fatherBT.block(0,0).data()->domainMapPtr()));
    else
        lagrangeMap.reset(new MAPEPETRA(*M_childBT.block(0,0).data()->domainMapPtr()));

    retVector.block(0).data().reset(new VECTOREPETRA(*lagrangeMap, LifeV::Unique));
    retVector.block(0).data()->zero();

    return retVector;
}

}
