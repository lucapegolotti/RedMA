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
setTimestep(double dt)
{
    AbstractAssembler::setTimestep(dt);
    // watch out: this might be a problem when using time adaptivity
    M_timeExtrapolator.setTimeStep(this->M_dt);
}

void
PseudoFSI::
setup()
{
    // first we need to setup the time extrapolator in order to have alpha for bdf
    // see boundary stiffness matrix
    unsigned int orderBDF = M_datafile("structure/time_integration_order", 1);
    M_timeExtrapolator.setBDForder(orderBDF);
    M_timeExtrapolator.setMaximumExtrapolationOrder(orderBDF);

    computeLameConstants();

    NavierStokesAssembler::setup();

    Vector velocityInitial(M_velocityFESpace->map());
    velocityInitial.zero();
    std::vector<Vector> initialStateVelocity;
    for (unsigned int i = 0; i < orderBDF; ++i)
        initialStateVelocity.push_back(velocityInitial);

    M_timeExtrapolator.initialize(initialStateVelocity);
    computeBoundaryIndicator();
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
    // VectorPtr aux(new Vector(M_velocityFESpace->map()));
    // aux->zero();
    M_boundaryIndicator.reset(new Vector(M_velocityFESpace->map()));
    M_boundaryIndicator->zero();
    // *aux += 1.0;
    // *M_boundaryIndicator = *M_boundaryMass * (*aux);
    // *M_boundaryIndicator = (*M_boundaryIndicator != 0);
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
                dt / M_timeExtrapolator.alpha() * (
                2  *  M_lameII *
                0.5 * dot(
                (grad(phi_j) - grad(phi_j)*outerProduct(Nface, Nface))
                + transpose(grad(phi_j) - grad(phi_j)*outerProduct(Nface, Nface)),
                (grad (phi_i) - grad(phi_i)*outerProduct(Nface, Nface))) +
                M_lameI *
                dot(value(Eye),(grad(phi_j) - grad(phi_j)*outerProduct(Nface, Nface))) *
                dot(value(Eye),(grad(phi_i) - grad(phi_i)*outerProduct(Nface, Nface))))
              ) >>  M_boundaryStiffness;

    M_boundaryStiffness->globalAssemble();
    *M_A += *M_boundaryStiffness;
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
    std::vector<VectorPtr> retVec = NavierStokesAssembler::computeF();

    VectorPtr rhsDisplacement(new Vector(M_velocityFESpace->map()));
    M_timeExtrapolator.rhsContribution(*rhsDisplacement);

    *retVec[0] -= *M_boundaryStiffness * (*rhsDisplacement);

    return retVec;
}

// attention: here we assume that the solutions in prevSolutions are the ones
// we have to use to update the displacement field..
void
PseudoFSI::
postProcess()
{
    VectorPtr rhsDisplacement(new Vector(M_velocityFESpace->map()));
    VectorPtr curDisplacement(new Vector(M_velocityFESpace->map()));

    M_timeExtrapolator.rhsContribution(*rhsDisplacement);

    curDisplacement->zero();
    *curDisplacement = (*M_prevSolution[0]) * (*M_boundaryIndicator);
    *curDisplacement += *rhsDisplacement;
    *curDisplacement *= (this->M_dt / M_timeExtrapolator.alpha());
    *M_displacementExporter = *curDisplacement;
    M_timeExtrapolator.shift(*curDisplacement);
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
