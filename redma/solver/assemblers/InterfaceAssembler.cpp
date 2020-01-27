#include "InterfaceAssembler.hpp"

namespace RedMA
{

template <>
SHP(LifeV::QuadratureRule)
InterfaceAssembler<VectorEp, MatrixEp>::
generateQuadratureRule(std::string tag) const
{
    using namespace LifeV;

    SHP(QuadratureRule) customRule;
    if (!std::strcmp(tag.c_str(),"STRANG10"))
    {
        // rule taken from https://people.sc.fsu.edu/~jburkardt/datasets/quadrature_rules_tri/quadrature_rules_tri.html
        customRule.reset(new QuadratureRule("STRANG10",TRIANGLE, 3, 7, 0));

        QuadraturePoint p1(0.333333333333333,0.333333333333333,0,
                           -0.149570044467670/2);
        customRule->addPoint(p1);

        QuadraturePoint p2(0.479308067841923,0.260345966079038,0,
                           0.175615257433204/2);
        customRule->addPoint(p2);

        QuadraturePoint p3(0.260345966079038,0.479308067841923,0,
                           0.175615257433204/2);
        customRule->addPoint(p3);

        QuadraturePoint p4(0.260345966079038,0.260345966079038,0,
                           0.175615257433204/2);
        customRule->addPoint(p4);

        QuadraturePoint p5(0.869739794195568,0.065130102902216,0,
                           0.053347235608839/2);
        customRule->addPoint(p5);

        QuadraturePoint p6(0.065130102902216,0.869739794195568,0,
                           0.053347235608839/2);
        customRule->addPoint(p6);

        QuadraturePoint p7(0.065130102902216,0.065130102902216,0,
                           0.053347235608839/2);
        customRule->addPoint(p7);

        QuadraturePoint p8(0.638444188569809,0.312865496004875,0,
                           0.077113760890257/2);
        customRule->addPoint(p8);

        QuadraturePoint p9(0.638444188569809,0.048690315425316,0,
                           0.077113760890257/2);
        customRule->addPoint(p9);

        QuadraturePoint p10(0.312865496004875,0.638444188569809,0,
                           0.077113760890257/2);
        customRule->addPoint(p10);

        QuadraturePoint p11(0.312865496004875,0.048690315425316,0,
                           0.077113760890257/2);
        customRule->addPoint(p11);

        QuadraturePoint p12(0.048690315425316,0.638444188569809,0,
                           0.077113760890257/2);
        customRule->addPoint(p12);

        QuadraturePoint p13(0.048690315425316,0.312865496004875,0,
                           0.077113760890257/2);
        customRule->addPoint(p13);
    }
    else
    {
        throw new Exception("Quadrature rule " + tag + " not implemented!");
    }
    return customRule;
}

template <>
void
InterfaceAssembler<VectorEp, MatrixEp>::
buildMapLagrange(SHP(BasisFunctionFunctor) bfs)
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
buildCouplingVectors(SHP(BasisFunctionFunctor) bfs,
                     const GeometricFace& face,
                     SHP(aAssembler<VectorEp COMMA MatrixEp>) assembler) const
{
    using namespace LifeV;
    using namespace ExpressionAssembly;

    QuadratureBoundary boundaryQuadRule(buildTetraBDQR
                                       (*generateQuadratureRule("STRANG10")));
    std::vector<VectorEp> couplingVectors;

    unsigned int nBasisFunctions = bfs->getNumBasisFunctions();

    SHP(ETFESPACE3) etfespace = assembler->getETFESpaceCoupling();
    MAPEPETRA map = etfespace->map();
    SHP(MESH) mesh = etfespace->mesh();
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
            SHP(VECTOREPETRA) currentMode(new VECTOREPETRA(map, LifeV::Repeated));

            bfs->setIndex(i);
            integrate(boundary(mesh, faceFlag),
                      boundaryQuadRule,
                      etfespace,
                      eval(bfs, X) * dot(versor, phi_i)
                  ) >> currentMode;
            couplingVectors[count].data() = currentMode;
            count++;
        }
    }

    return couplingVectors;
}

template <>
void
InterfaceAssembler<VectorEp, MatrixEp>::
buildCouplingMatrices(SHP(AssemblerType) assembler, const GeometricFace& face,
                      BlockMatrix<MatrixEp>& matrixT, BlockMatrix<MatrixEp>& matrix)
{
    SHP(BasisFunctionFunctor) bfs;
    bfs =  BasisFunctionFactory(M_data.getDatafile(), face);

    std::vector<VectorEp> couplingVectors;
    couplingVectors = buildCouplingVectors(bfs, face, assembler);

    matrixT.resize(assembler->getNumComponents(), 1);
    matrix.resize(1, assembler->getNumComponents());
    // we assume that the first field is the one to be coupled
    matrixT.block(0,0).softCopy(MatrixEp(couplingVectors));
    matrix.block(0,0).softCopy(matrixT.block(0,0).transpose());

    assembler->getBCManager()->apply0DirichletMatrix(matrixT,
                                                     assembler->getFESpaceBCs(),
                                                     assembler->getComponentBCs(),
                                                     0.0);
}

template <>
void
InterfaceAssembler<VectorEp, MatrixEp>::
buildCouplingMatrices()
{
    unsigned int indexOutlet = M_interface.M_indexOutlet;

    auto asFather = M_interface.M_assemblerFather;
    GeometricFace outlet = asFather->getTreeNode()->M_block->getOutlet(indexOutlet);

    buildCouplingMatrices(asFather, outlet, M_fatherBT, M_fatherB);

    auto asChild = M_interface.M_assemblerChild;
    GeometricFace inlet = asChild->getTreeNode()->M_block->getInlet();
    // I invert the normal of the face such that it is the same as the outlet
    inlet.M_normal *= (-1.);

    buildCouplingMatrices(asChild, inlet, M_childBT, M_childB);
}
}
