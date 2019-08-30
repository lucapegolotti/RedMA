#include <PseudoFSI.hpp>

namespace RedMA
{

PseudoFSI::
PseudoFSI(const GetPot& datafile, commPtr_Type comm,
          const TreeNodePtr& treeNode, bool verbose) :
  NavierStokesAssembler(datafile, comm, treeNode, verbose)
{
    // this allows us to do not set anything at the wall
    M_addNoslipBC = false;
}

void
PseudoFSI::
setup()
{
    computeLameConstants();

    NavierStokesAssembler::setup();

    computeBoundaryIndicator();

    // set mass displacement to identity
    M_massDisplacement.reset(new Matrix(M_velocityFESpace->map()));
    M_massDisplacement->insertOneDiagonal();
    M_massDisplacement->globalAssemble();

    // adding map for displacement
    M_primalMaps.push_back(M_velocityFESpace->mapPtr());
}

void
PseudoFSI::
assembleMassMatrix()
{
    using namespace LifeV::ExpressionAssembly;

    NavierStokesAssembler::assembleMassMatrix();

    printlog(YELLOW, "Assembling boundary mass matrix ...\n", M_verbose);

    LifeV::QuadratureBoundary myBDQR(LifeV::buildTetraBDQR(LifeV::quadRuleTria4pt));

    M_boundaryMass.reset(new Matrix(M_velocityFESpace->map()));
    M_boundaryMass->zero();

    double density = M_datafile("structure/density", 1.2);
    double thickness = M_datafile("structure/thickness", 0.1);
    unsigned int wallFlag = M_datafile("structure/flag", 10);
    integrate(boundary(M_velocityFESpaceETA->mesh(), wallFlag),
               myBDQR,
               M_velocityFESpaceETA,
               M_velocityFESpaceETA,
               value(density * thickness) * dot (phi_i, phi_j)
              ) >> M_boundaryMass;

    M_boundaryMass->globalAssemble();
    *M_M += *M_boundaryMass;
}

void
PseudoFSI::
computeBoundaryIndicator()
{
    M_boundaryIndicator.reset(new Vector(M_velocityFESpace->map()));
    M_boundaryIndicator->zero();

    LifeV::BCFunctionBase oneFunction(fOne);

    BoundaryConditionPtr bcs;
    bcs.reset(new LifeV::BCHandler);

    const unsigned int wallFlag = 10;

    bcs->addBC("Wall", wallFlag, LifeV::Essential,
                LifeV::Full, oneFunction, 3);

    updateBCs(bcs, M_velocityFESpace);

    bcManageRhs(*M_boundaryIndicator, *M_velocityFESpace->mesh(), M_velocityFESpace->dof(),
                *bcs, M_velocityFESpace->feBd(), 1.0, 0.0);
}

NavierStokesAssembler::MatrixPtr
PseudoFSI::
getMassMatrix(const unsigned int& blockrow,
              const unsigned int& blockcol)
{
    if (blockrow == 0 && blockcol == 0)
    {
        // we copy the matrix so that we cannot change M_M (e.g. by applying
        // boundary conditions)
        MatrixPtr returnMatrix(new Matrix(*M_M));
        return returnMatrix;
    }
    else if (blockrow == 2 && blockcol == 2)
    {
        MatrixPtr returnMatrix(new Matrix(*M_massDisplacement));
        return returnMatrix;
    }

    return nullptr;
}


NavierStokesAssembler::MatrixPtr
PseudoFSI::
getJacobian(const unsigned int& blockrow, const unsigned int& blockcol)
{
    if (blockrow < 2 && blockcol < 2)
        return NavierStokesAssembler::getJacobian(blockrow,blockcol);
    if (blockrow == 2 && blockcol == 0)
    {
        MatrixPtr returnMatrix(new Matrix(*M_massDisplacement));
        // *returnMatrix *= (-1.0);
        returnMatrix->globalAssemble();
        return returnMatrix;
    }
    if (blockrow == 0 && blockcol == 2)
    {
        MatrixPtr returnMatrix(new Matrix(*M_boundaryStiffness));
        *returnMatrix *= (-1.0);
        returnMatrix->globalAssemble();
        return returnMatrix;
    }
    return nullptr;
}

void
PseudoFSI::
assembleStiffnessMatrix()
{
    using namespace LifeV::ExpressionAssembly;

    NavierStokesAssembler::assembleStiffnessMatrix();

    printlog(YELLOW, "Assembling boundary stiffness matrix ...\n", M_verbose);

    LifeV::QuadratureBoundary myBDQR(LifeV::buildTetraBDQR(LifeV::quadRuleTria4pt));

    M_boundaryStiffness.reset(new Matrix(M_velocityFESpace->map()));
    M_boundaryStiffness->zero();
    LifeV::MatrixSmall<3, 3> Eye;
    Eye *= 0.0;
    Eye[0][0] = 1;
    Eye[1][1] = 1;
    Eye[2][2] = 1;
    unsigned int wallFlag = M_datafile("structure/flag", 10);
    // not so nice: we rely on the datafile because dt has not been set yet
    // in the extrapolator. I might think of another solution if using
    // adapativity in time
    double  dt = M_datafile("time_discretization/dt", 0.01);

    integrate(boundary(M_velocityFESpaceETA->mesh(), wallFlag),
                myBDQR,
                M_velocityFESpaceETA,
                M_velocityFESpaceETA,
                ( 2  *  M_lameII *
                0.5 * dot(
                (grad(phi_j) - grad(phi_j)*outerProduct(Nface, Nface))
                + transpose(grad(phi_j) - grad(phi_j)*outerProduct(Nface, Nface)),
                (grad (phi_i) - grad(phi_i)*outerProduct(Nface, Nface))) +
                M_lameI *
                dot(value(Eye),(grad(phi_j) - grad(phi_j)*outerProduct(Nface, Nface))) *
                dot(value(Eye),(grad(phi_i) - grad(phi_i)*outerProduct(Nface, Nface))))
              ) >>  M_boundaryStiffness;

    M_boundaryStiffness->globalAssemble();
}

void
PseudoFSI::
computeLameConstants()
{
    double poisson = M_datafile("structure/poisson", 0.45);
    double young = M_datafile("structure/young", 4e6);
    double thickness = M_datafile("structure/thickness", 0.1);

    M_lameI = (thickness * young * poisson)/((1. - 2*poisson)*(1. + poisson));
    M_lameII = thickness * young/(2. * (1. + poisson));
}

std::vector<PseudoFSI::VectorPtr>
PseudoFSI::
computeF()
{
    std::vector<VectorPtr> retVec;
    std::vector<VectorPtr> retVecNS = NavierStokesAssembler::computeF();

    // N.B: we want the coupling part to be in the last position of the residual
    *retVecNS[0] -= *M_boundaryStiffness * (*M_prevSolution[2]);

    retVec.push_back(retVecNS[0]);
    retVec.push_back(retVecNS[1]);

    VectorPtr F3(new Vector(M_velocityFESpace->map()));
    *F3 += *M_prevSolution[0];
    retVec.push_back(F3);
    for (unsigned int i = 2; i < retVecNS.size(); i++)
        retVec.push_back(retVecNS[i]);

    return retVec;
}

std::vector<NavierStokesAssembler::VectorPtr>
PseudoFSI::
computeFder()
{
    std::vector<VectorPtr> retVec;
    std::vector<VectorPtr> retVecNS = NavierStokesAssembler::computeF();

    retVec.push_back(retVecNS[0]);
    retVec.push_back(retVecNS[1]);

    VectorPtr F3(new Vector(M_velocityFESpace->map()));
    F3->zero();

    retVec.push_back(F3);
    for (unsigned int i = 2; i < retVecNS.size(); i++)
        retVec.push_back(retVecNS[i]);

    return retVec;
}

// attention: here we assume that the solutions in prevSolutions are the ones
// we have to use to update the displacement field..
void
PseudoFSI::
postProcess()
{
    *M_displacementExporter = *M_prevSolution[2] * (*M_boundaryIndicator);
}

void
PseudoFSI::
setExporter()
{
    NavierStokesAssembler::setExporter();
    M_displacementExporter.reset(new Vector(M_velocityFESpace->map(),
                                            M_exporter->mapType()));

    M_exporter->addVariable(LifeV::ExporterData<Mesh>::VectorField,
                         "displacement", M_velocityFESpace, M_displacementExporter, 0.0);
}

}  // namespace RedMA
