// implementation of template class

namespace RedMA
{

template <class AssemblerType>
TimeMarchingAlgorithm<AssemblerType>::
TimeMarchingAlgorithm(const GetPot& datafile,
                      GlobalAssemblerType* assembler,
                      commPtr_Type comm,
                      bool verbose) :
  M_datafile(datafile),
  M_globalAssembler(assembler),
  M_comm(comm),
  M_verbose(verbose)
{
    M_solution.reset(new Vector(assembler->getGlobalMap()));
    M_solution->zero();
}

template <class AssemblerType>
void
TimeMarchingAlgorithm<AssemblerType>::
solveLinearSystem(MatrixPtr matrix, VectorPtr rhs, VectorPtr sol)
{
    // reset solution vector
    sol->zero();

    // solver part
    LifeV::LinearSolver linearSolver(M_comm);
    linearSolver.setOperator(matrix);

    Teuchos::RCP<Teuchos::ParameterList> aztecList =
                                       Teuchos::rcp(new Teuchos::ParameterList);

    std::string xmlSolverData = M_datafile("solver/XMLdatafile",
                                           "SolverParamList.xml");

    aztecList = Teuchos::getParametersFromXmlFile(xmlSolverData);
    linearSolver.setParameters(*aztecList);

    typedef LifeV::PreconditionerIfpack             precIfpack;
    // typedef std::shared_ptr<precIfpack>        precIfpackPtr;
    precIfpack* precRawPtr = new precIfpack;

    // we set to look for the "fake" precMLL entry in order to set the
    // default parameters of ML preconditioner
    GetPot dummyDatafile;
    precRawPtr->setDataFromGetPot(dummyDatafile, "precMLL");
    std::shared_ptr<LifeV::Preconditioner> precPtr;
    precPtr.reset(precRawPtr);
    precPtr->parametersList().set( "precType", "Amesos" );
    precPtr->parametersList().set( "overlap level", 2 );
    precPtr->parametersList().set( "amesos: solver type",
             M_datafile( "prec/ifpack/amesos/solvertype", "Amesos_Umfpack" ) );
    // print information about the solver
    // precPtr->parametersList().print ( std::cout );

    linearSolver.setPreconditioner(precPtr);
    linearSolver.setRightHandSide(rhs);
    linearSolver.solve(sol);
}

template <class AssemblerType>
typename TimeMarchingAlgorithm<AssemblerType>::VectorPtr
TimeMarchingAlgorithm<AssemblerType>::
getSolution()
{
    return M_solution;
}

}  // namespace RedMA
