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

    setSolversOptions();
}

template <class AssemblerType>
void
TimeMarchingAlgorithm<AssemblerType>::
solveLinearSystem(GlobalBlockMatrix matrix,
                  VectorPtr rhs, VectorPtr sol,
                  GlobalBlockMatrix* matrixPrec)
{
    // reset solution vector
    sol->zero();

    auto grid = matrix.getGrid();
    M_oper.reset(new LifeV::Operators::GlobalSolverOperator());
    M_oper->setUp(grid, M_comm);
    if (matrixPrec == nullptr)
        buildPreconditioner(matrix);
    else
        buildPreconditioner(*matrixPrec);

    std::string solverType(M_pListLinSolver->get<std::string>("Linear Solver Type"));
    M_invOper.reset(LifeV::Operators::
                InvertibleOperatorFactory::instance().createObject(solverType));

    M_invOper->setParameterList(M_pListLinSolver->sublist(solverType));
    M_invOper->setOperator(M_oper);
    M_invOper->setPreconditioner(M_prec);

    M_invOper->ApplyInverse(rhs->epetraVector(), sol->epetraVector());
}

template <class AssemblerType>
void
TimeMarchingAlgorithm<AssemblerType>::
buildPreconditioner(GlobalBlockMatrix matrix)
{
    AssemblerType as(M_datafile, M_comm, nullptr, false);
    if (as.numberOfBlocks() == 2)
    {
        bool steady = M_datafile("time_discretization/steady", false);
        if (steady)
            M_prec.reset(new LifeV::Operators::GlobalSIMPLEOperator("STEADY"));
        else
            M_prec.reset(new LifeV::Operators::GlobalSIMPLEOperator("SIMPLE"));
    }
    else
        M_prec.reset(new LifeV::Operators::GlobalSIMPLEOperatorPseudoFSI());
    M_prec->setVerbose(M_verbose);
    M_prec->setSolversOptions(*M_solversOptions);
    M_prec->setUp(matrix, M_comm);
}

template <class AssemblerType>
unsigned int
TimeMarchingAlgorithm<AssemblerType>::
getOrder()
{
    return M_order;
}

template <class AssemblerType>
typename TimeMarchingAlgorithm<AssemblerType>::VectorPtr
TimeMarchingAlgorithm<AssemblerType>::
getSolution()
{
    return M_solution;
}

template <class AssemblerType>
void
TimeMarchingAlgorithm<AssemblerType>::
setInitialCondition(VectorPtr initalCondition)
{
    M_solution = initalCondition;
}

template <class AssemblerType>
void
TimeMarchingAlgorithm<AssemblerType>::
setSolversOptions()
{
    std::string optionsPrec = M_datafile("fluid/options_preconditioner",
                                         "solversOptionsFast");
    optionsPrec += ".xml";
    M_solversOptions = Teuchos::getParametersFromXmlFile(optionsPrec);
    std::shared_ptr<Teuchos::ParameterList> monolithicOptions;
    monolithicOptions.reset(
        new Teuchos::ParameterList(M_solversOptions->sublist("MonolithicOperator")));
    M_pListLinSolver = monolithicOptions;
}

template <class AssemblerType>
void
TimeMarchingAlgorithm<AssemblerType>::
solveNonLinearSystem(std::function<VectorPtr(VectorPtr)> fun,
                     std::function<GlobalBlockMatrix(VectorPtr)> jac,
                     VectorPtr sol,
                     const double& tol, const unsigned int& itMax,
                     std::function<GlobalBlockMatrix(VectorPtr)>* jacPrec)
{
    VectorPtr incr(new Vector(sol->map()));
    double err = tol + 1;
    unsigned int count = 1;
    VectorPtr curF;
    while (err > tol && count < itMax)
    {
        // VectorPtr curF(new Vector(x0->map()));
        curF = fun(sol);
        err = curF->norm2();
        std::string msg("[SolveNonLinearSystem]");
        msg += " solving, iteration = " + std::to_string(count) + ", ";
        std::ostringstream streamOb;
        streamOb << err;
        msg += " error = " + streamOb.str() + "\n";
        printlog(YELLOW, msg, M_verbose);
        if (err > tol)
        {
            incr->zero();
            GlobalBlockMatrix curJac = jac(sol);
            if (jacPrec == nullptr)
                solveLinearSystem(curJac, curF, incr);
            else
            {
                std::function<GlobalBlockMatrix(VectorPtr)> jacPrecRaw = *jacPrec;
                GlobalBlockMatrix curJacPrec = jacPrecRaw(sol);
                curJacPrec.printPattern();
                solveLinearSystem(curJac, curF, incr, &curJacPrec);
            }
            *sol -= *incr;
            count++;
        }
    }
    // *sol = *curF;
    if (count != itMax)
    {
        std::string msg("[SolveNonLinearSystem]");
        msg += " convergence, iteration = " + std::to_string(count) + ", ";
        msg += " error = " + std::to_string(err) + "\n";
        printlog(YELLOW, msg, M_verbose);
    }
    else
    {
        std::string msg("[SolveNonLinearSystem]");
        msg += " did not reach convergence!\n";
        printlog(RED, msg, M_verbose);
    }
}

}  // namespace RedMA
