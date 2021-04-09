#include "NonAffineDeformer.hpp"

namespace RedMA
{

NonAffineDeformer::
NonAffineDeformer(shp<MESH> mesh, EPETRACOMM comm, bool verbose) :
  M_mesh(mesh),
  M_comm(comm),
  M_verbose(verbose)
{
    // hardcoded values for the moment
    const double young = 100;// 3e6;
    const double poisson = 0.3;

    M_fespace.reset(new FESPACE(M_mesh, "P1", 3, M_comm));
	M_fespaceETA.reset(new ETFESPACE3(M_fespace->mesh(),
                                     &(M_fespace->refFE()), M_comm));

    M_stiffness.reset(new MATRIXEPETRA(M_fespace->map()));
    assembleStiffness(young, poisson);

    M_rhs.reset(new VECTOREPETRA(M_fespace->map(), LifeV::Unique));
    M_rhs->zero();
}

void
NonAffineDeformer::
assembleStiffness(const double young, const double poisson)
{
    using namespace LifeV::ExpressionAssembly;

    // compute lame coefficients
    double lambda = (young * poisson)/((1.0 + poisson) * (1.0 - 2.0 * poisson));
    double mu = young / (2.0 * (1.0 + poisson));

    integrate ( elements (M_fespaceETA->mesh() ),
                M_fespace->qr(), M_fespaceETA, M_fespaceETA,
                value(mu) * dot(grad(phi_i), grad(phi_j) + transpose(grad(phi_j))) +
                value(lambda) * div(phi_i) * div(phi_j)
            ) >> M_stiffness;

    M_stiffness->globalAssemble();
}

void
NonAffineDeformer::
applyBCs(shp<LifeV::BCHandler> bcs)
{
    bcs->bcUpdate(*M_fespace->mesh(), M_fespace->feBd(), M_fespace->dof());

    LifeV::bcManage(*M_stiffness, *M_rhs, *M_fespace->mesh(), M_fespace->dof(),
                    *bcs, M_fespace->feBd(), 1.0);
}

void
NonAffineDeformer::
deformMesh(LifeV::MeshUtility::MeshTransformer<MESH>& transformer)
{
    shp<VECTOREPETRA> displacement = solveSystem();
    transformer.moveMesh(*displacement,  M_fespace->dof().numTotalDof());
}

shp<VECTOREPETRA>
NonAffineDeformer::
solveSystem()
{
    shp<VECTOREPETRA> solution;
    // reset solution vector
    solution.reset(new VECTOREPETRA(M_fespace->map(), LifeV::Unique));
    solution->zero();

    // solver part
    LifeV::LinearSolver linearSolver(M_comm);
    linearSolver.setOperator(M_stiffness);

    Teuchos::RCP<Teuchos::ParameterList> aztecList =
                                       Teuchos::rcp(new Teuchos::ParameterList);
    aztecList = Teuchos::getParametersFromXmlFile(M_XMLsolver);
    linearSolver.setParameters(*aztecList);

    typedef LifeV::PreconditionerML         precML_type;
    typedef shp<precML_type>                precMLPtr_type;
    precML_type * precRawPtr;
    precRawPtr = new precML_type;
    // we set to look for the "fake" precMLL entry in order to set the
    // default parameters of ML preconditioner
    GetPot dummyDatafile;
    precRawPtr->setDataFromGetPot(dummyDatafile, "precMLL");
    shp<LifeV::Preconditioner> precPtr;
    precPtr.reset(precRawPtr);

    linearSolver.setPreconditioner(precPtr);
    linearSolver.setRightHandSide(M_rhs);
    linearSolver.solve(solution);

    return solution;
}

void
NonAffineDeformer::
setXMLsolver(std::string filename)
{
    M_XMLsolver = filename;
}

}  // namespace RedMA
