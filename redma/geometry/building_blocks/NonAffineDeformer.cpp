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
    const double young = 100;
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
    shp<VECTOREPETRA> displacement = solveSystem("Ifpack");
    transformer.moveMesh(*displacement, M_fespace->dof().numTotalDof());
}

void
NonAffineDeformer::
deformMeshComposite(LifeV::MeshUtility::MeshTransformer<MESH>& transformer, shp<VECTOREPETRA> displacement)
{
    transformer.moveMesh(*displacement, M_fespace->dof().numTotalDof());
}

shp<VECTOREPETRA>
NonAffineDeformer::
solveSystem(const std::string& precType)
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

//    std::string optionsPrec = M_data("preconditioner/options",
//                                     "datafiles/solversOptionsFast");
//    optionsPrec += ".xml";
//    Teuchos::RCP<Teuchos::ParameterList> precList = Teuchos::getParametersFromXmlFile(optionsPrec);
//    shp<Teuchos::ParameterList> precOptions;
//    precOptions.reset(
//            new Teuchos::ParameterList(M_solversOptionsInner->sublist(precType.c_str())));
//    linearSolver.setParameters(*precOptions);

    shp<LifeV::Preconditioner> precPtr;

    if (!std::strcmp(precType.c_str(), "ML")) {
        typedef LifeV::PreconditionerML precML_type;
        typedef shp<precML_type>  precMLPtr_type;
        precML_type* precRawPtr;
        precRawPtr = new precML_type;
        GetPot dummyDatafile;
        precRawPtr->setDataFromGetPot(dummyDatafile, "prec");
        precPtr.reset(precRawPtr);
    }
    else if (!std::strcmp(precType.c_str(), "Ifpack")) {
        typedef LifeV::PreconditionerIfpack precIf_type;
        typedef shp<precIf_type>  precIfPtr_type;
        precIf_type* precRawPtr;
        precRawPtr = new precIf_type;
        GetPot dummyDatafile;
        precRawPtr->setDataFromGetPot(dummyDatafile, "prec");
        precPtr.reset(precRawPtr);
    }
    else {
        throw new Exception("Unrecognized preconditioner type " + precType);
    }

    // linearSolver.setPreconditionerFromGetPot("datafiles/data", "preconditioner/deformation");
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
