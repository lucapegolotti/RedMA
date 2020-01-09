#include <TimeMarchingAlgorithm.hpp>

namespace RedMA
{

TimeMarchingAlgorithm::
TimeMarchingAlgorithm(const GetPot& datafile,
                      AbstractAssembler* assembler,
                      commPtr_Type comm,
                      bool verbose) :
  M_datafile(datafile),
  M_assembler(assembler),
  M_comm(comm),
  M_verbose(verbose)
{
    // M_solution.reset(new Vector(assembler->getGlobalMap()));
    // M_solution->zero();
    //
    // setSolversOptions();
}

void
TimeMarchingAlgorithm::
solveLinearSystem(BlockMatrix matrix,
                  VectorPtr rhs, VectorPtr sol,
                  BlockMatrix* matrixPrec)
{
    // // reset solution vector
    // sol->zero();
    //
    // auto grid = matrix.getGrid();
    // M_oper.reset(new LifeV::Operators::GlobalSolverOperator());
    // M_oper->setUp(grid, M_comm);
    // if (matrixPrec == nullptr)
    //     buildPreconditioner(matrix);
    // else
    //     buildPreconditioner(*matrixPrec);
    //
    // std::string solverType(M_pListLinSolver->get<std::string>("Linear Solver Type"));
    // M_invOper.reset(LifeV::Operators::
    //             InvertibleOperatorFactory::instance().createObject(solverType));
    //
    // M_invOper->setParameterList(M_pListLinSolver->sublist(solverType));
    // M_invOper->setOperator(M_oper);
    // M_invOper->setPreconditioner(M_prec);
    //
    // M_invOper->ApplyInverse(rhs->epetraVector(), sol->epetraVector());
}

void
TimeMarchingAlgorithm::
buildPreconditioner(BlockMatrix matrix)
{
    // std::shared_ptr<AbstractAssembler> as = AssemblersFactory(M_datafile, M_comm,
    //                                                           nullptr, M_verbose);
    // if (as->numberOfBlocks() == 2)
    // {
    //     bool steady = M_datafile("time_discretization/steady", false);
    //     bool exactSolveBlocks = M_datafile("linear_solve/exact_solve_primal", false);
    //     // if -1, all iterations are exact
    //     int numIterationsExactSolve = M_datafile("linear_solve/num_exact_iterations", 2);
    //
    //     if (numIterationsExactSolve > 0 && exactSolveBlocks)
    //     {
    //       throw new Exception("This feature is buggy!");
    //     }
    //     // we don't allow for exact solve in steady case. The problem is that we don't pass
    //     // the matrix to be solved to the preconditioner
    //     if (steady)
    //         M_prec.reset(new LifeV::Operators::GlobalSIMPLEOperator("STEADY", false, -1));
    //     else
    //         M_prec.reset(new LifeV::Operators::GlobalSIMPLEOperator("SIMPLE", exactSolveBlocks,
    //                                                                 numIterationsExactSolve));
    // }
    // else
    //     M_prec.reset(new LifeV::Operators::GlobalSIMPLEOperatorPseudoFSI());
    // M_prec->setVerbose(M_verbose);
    // M_prec->setSolversOptions(*M_solversOptions);
    // M_prec->setUp(matrix, M_comm);
}

unsigned int
TimeMarchingAlgorithm::
getOrder()
{
    return M_order;
}

BlockVector
TimeMarchingAlgorithm::
getSolution()
{
    return M_solution;
}

void
TimeMarchingAlgorithm::
setInitialCondition(BlockVector initalCondition)
{
    M_solution = initalCondition;
}

void
TimeMarchingAlgorithm::
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

void
TimeMarchingAlgorithm::
solveNonLinearSystem(std::function<VectorPtr(VectorPtr)> fun,
                     std::function<BlockMatrix(VectorPtr)> jac,
                     VectorPtr sol,
                     const double& tol, const unsigned int& itMax,
                     std::function<BlockMatrix(VectorPtr)>* jacPrec)
{
    // VectorPtr incr(new Vector(sol->map()));
    // double err = tol + 1;
    // unsigned int count = 1;
    // VectorPtr curF;
    // while (err > tol && count < itMax)
    // {
    //     // VectorPtr curF(new Vector(x0->map()));
    //     curF = fun(sol);
    //     err = curF->norm2();
    //     std::string msg("[SolveNonLinearSystem]");
    //     msg += " solving, iteration = " + std::to_string(count) + ", ";
    //     std::ostringstream streamOb;
    //     streamOb << err;
    //     msg += " error = " + streamOb.str() + "\n";
    //     printlog(YELLOW, msg, M_verbose);
    //     if (err > tol)
    //     {
    //         incr->zero();
    //         BlockMatrix curJac = jac(sol);
    //         if (jacPrec == nullptr)
    //             solveLinearSystem(curJac, curF, incr);
    //         else
    //         {
    //             std::function<BlockMatrix(VectorPtr)> jacPrecRaw = *jacPrec;
    //             BlockMatrix curJacPrec = jacPrecRaw(sol);
    //             solveLinearSystem(curJac, curF, incr, &curJacPrec);
    //         }
    //         *sol -= *incr;
    //         count++;
    //     }
    // }
    // // *sol = *curF;
    // if (count != itMax)
    // {
    //     std::string msg("[SolveNonLinearSystem]");
    //     msg += " convergence, iteration = " + std::to_string(count) + ", ";
    //     msg += " error = " + std::to_string(err) + "\n";
    //     printlog(YELLOW, msg, M_verbose);
    // }
    // else
    // {
    //     std::string msg("[SolveNonLinearSystem]");
    //     msg += " did not reach convergence!\n";
    //     printlog(RED, msg, M_verbose);
    // }
}

}  // namespace RedMA
