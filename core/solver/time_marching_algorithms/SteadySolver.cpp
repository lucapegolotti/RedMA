#include <SteadySolver.hpp>

namespace RedMA
{

SteadySolver::
SteadySolver(const GetPot& datafile,
             AbstractAssembler* assembler,
             commPtr_Type comm,
             bool verbose) :
  TimeMarchingAlgorithm(datafile, assembler, comm, verbose)
{
    this->M_order = 0;
}

void
SteadySolver::
solveTimestep(const double &time, double &dt)
{
    // std::string msg("[SteadySolver] solving ... \n");
    // printlog(MAGENTA, msg, M_verbose);
    // M_globalAssembler->applyBCsVector(M_solution, 1.0, time + dt,
    //                                   &AbstractAssembler::applyBCsBackwardEuler);
    // std::function<VectorPtr(VectorPtr)> mFun =
    //       std::bind(&SteadySolver::assembleF,
    //                 this, time + dt, std::placeholders::_1, dt);
    //
    // std::function<BlockMatrix(VectorPtr)> mJac =
    //       std::bind(&SteadySolver::assembleJac,
    //                 this, time + dt, std::placeholders::_1, dt);
    //
    // std::function<BlockMatrix(VectorPtr)> mJacPrec =
    //       std::bind(&SteadySolver::assembleJacPrec,
    //                 this, time + dt, std::placeholders::_1, dt);
    //
    // double tol = M_datafile("newton_method/tol", 1e-5);
    // double maxit = M_datafile("newton_method/maxit", 5);
    // this->solveNonLinearSystem(mFun, mJac, M_solution, tol, maxit, &mJacPrec);
    //
    // printlog(MAGENTA, "done\n", M_verbose);
}

typename SteadySolver::VectorPtr
SteadySolver::
assembleF(const double& time, BlockVector tentativeSol, const double& dt)
{
    // std::string msg("assembling F\n");
    // printlog(GREEN, msg, M_verbose);
    // M_assembler->setTimeAndPrevSolution(time, tentativeSol);
    // VectorPtr retF(new Vector(M_solution->map()));
    // retF->zero();
    //
    // *retF += *(M_assembler->computeF_());
    // M_assembler->applyBCsVector(retF, 0.0, time,
    //                             &AbstractAssembler::applyBCsBackwardEuler);
    // return retF;
}

BlockMatrix
SteadySolver::
assembleJac(const double& time, BlockVector tentativeSol, const double& dt)
{
    // std::string msg("assembling Jacobian\n");
    // printlog(GREEN, msg, M_verbose);
    //
    // M_assembler->setTimeAndPrevSolution(time, tentativeSol, false);
    // double diagonalCoefficient = 1.0;
    // BlockMatrix retJac = M_assembler->getJacobianF(true, &diagonalCoefficient);
    //
    // return retJac;
}

BlockMatrix
SteadySolver::
assembleJacPrec(const double& time, BlockVector tentativeSol, const double& dt)
{
    // std::string msg("assembling preconditioner Jacobian\n");
    // printlog(GREEN, msg, M_verbose);
    // // we don't assemble the blocks again because this function is called
    // // after assembleF (and the stabilization blocks are already assembled there)
    // M_assembler->setTimeAndPrevSolution(time, tentativeSol, false);
    // double diagonalCoefficient = 1.0;
    // BlockMatrix retJac = M_assembler->getJacobianFprec(true, &diagonalCoefficient);
    // return retJac;
}

}  // namespace RedMA
