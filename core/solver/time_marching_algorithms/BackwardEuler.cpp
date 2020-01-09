#include <BackwardEuler.hpp>

namespace RedMA
{

BackwardEuler::
BackwardEuler(const GetPot& datafile,
              AbstractAssembler* assembler,
              commPtr_Type comm,
              bool verbose) :
  TimeMarchingAlgorithm(datafile, assembler, comm, verbose)
{
    double diagonalCoefficient = 1.0;
    // first we store the mass with no boundary conditions
    assembler->assembleMassMatrix(&diagonalCoefficient);

    M_prevSolution.reset(new Vector(assembler->getGlobalMap()));
    M_prevSolution->zero();

    this->M_order = 1;
}

void
BackwardEuler::
solveTimestep(const double &time, double &dt)
{
    // std::string msg("[BackwardEuler] solving, time = ");
    // msg += std::to_string(time) + " ...\n";
    // printlog(MAGENTA, msg, M_verbose);
    // // attention, here we are requesting a hard copy. So BlockVector must be provided
    // // with a hard copy member function
    // *M_prevSolution = *M_solution;
    // M_assembler->applyBCsBackwardEuler(M_solution, 1.0, time + dt);
    //
    // std::function<VectorPtr(VectorPtr)> mFun =
    //       std::bind(&BackwardEuler::assembleF,
    //                 this, time + dt, std::placeholders::_1, dt);
    //
    // std::function<BlockMatrix(VectorPtr)> mJac =
    //       std::bind(&BackwardEuler::assembleJac,
    //                 this, time + dt, std::placeholders::_1, dt);
    //
    // double tol = M_datafile("newton_method/tol", 1e-5);
    // double maxit = M_datafile("newton_method/maxit", 5);
    // this->solveNonLinearSystem(mFun, mJac, M_solution, tol, maxit);
    //
    // // M_assembler->checkResidual(M_solution, M_prevSolution, dt);
    //
    // printlog(MAGENTA, "done\n", M_verbose);
}

typename BackwardEuler::VectorPtr
BackwardEuler::
assembleF(const double& time, BlockVector tentativeSol, const double& dt)
{
    // std::string msg("assembling F\n");
    // printlog(GREEN, msg, M_verbose);
    // M_assembler->setTimeAndPrevSolution(time, tentativeSol);
    // VectorPtr retF(new Vector(M_solution->map()));
    // retF->zero();
    //
    // *retF = (M_assembler->getGlobalMass()) * (*M_prevSolution);
    // *retF *= (-1.0);
    // *retF += (M_assembler->getGlobalMass()) * (*tentativeSol);
    //
    // *retF -= (dt) * (*(M_assembler->computeF_()));
    // M_assembler->applyBCsVector(retF, 0.0, time,
    //                                   &AbstractAssembler::applyBCsBackwardEuler);
    //
    // return retF;
}

BlockMatrix
BackwardEuler::
assembleJac(const double& time, BlockVector tentativeSol, const double& dt)
{
    // std::string msg("assembling Jacobian\n");
    // printlog(GREEN, msg, M_verbose);
    // // we don't assemble the blocks again because this function is called
    // // after assembleF (and the stabilization blocks are already assembled there)
    // M_assembler->setTimeAndPrevSolution(time, tentativeSol, false);
    // double diagonalCoefficient = 0;
    // BlockMatrix retJac = M_assembler->getJacobianF(true,
    //                                                       &diagonalCoefficient);
    //
    // retJac *= (-dt);
    //
    // // trick to compute the jacobian of the part corresponding to H(un+1) * un
    // M_assembler->setTimeAndPrevSolution(time, M_prevSolution, false);
    // BlockMatrix massUn = M_assembler->getGlobalMassJacVelocity();
    //
    // massUn *= (-1);
    //
    // retJac.add(massUn);
    // retJac.add(M_assembler->getGlobalMassJac());
    // return retJac;
}

}  // namespace RedMA
