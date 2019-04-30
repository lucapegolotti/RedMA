// implementation of template class

namespace RedMA
{

template <class AssemblerType>
TimeMarchingAlgorithm<AssemblerType>::
TimeMarchingAlgorithm(const GetPot& datafile,
                      GlobalAssemblerType* assembler,
                      commPtr_Type comm) :
  M_datafile(datafile),
  M_globalAssembler(assembler),
  M_comm(comm)
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

    typedef LifeV::PreconditionerML         precML_type;
    typedef std::shared_ptr<precML_type>    precMLPtr_type;
    precML_type * precRawPtr;
    precRawPtr = new precML_type;
    // we set to look for the "fake" precMLL entry in order to set the
    // default parameters of ML preconditioner
    GetPot dummyDatafile;
    precRawPtr->setDataFromGetPot(dummyDatafile, "precMLL");
    std::shared_ptr<LifeV::Preconditioner> precPtr;
    precPtr.reset(precRawPtr);

    linearSolver.setPreconditioner(precPtr);
    linearSolver.setRightHandSide(rhs);
    linearSolver.solve(sol);
}

}  // namespace RedMA
