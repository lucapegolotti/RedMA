// implementation of template class

namespace RedMA
{

template <class AssemblerType>
SteadySolver<AssemblerType>::
SteadySolver(const GetPot& datafile,
             GlobalAssemblerType* assembler,
             commPtr_Type comm,
             bool verbose) :
  TimeMarchingAlgorithm<AssemblerType>(datafile, assembler, comm, verbose)
{
    this->M_order = 0;
}

template <class AssemblerType>
void
SteadySolver<AssemblerType>::
solveTimestep(const double &time, double &dt)
{
    std::string msg("[SteadySolver] solving ... \n");
    printlog(MAGENTA, msg, M_verbose);
    M_globalAssembler->applyBCsVector(M_solution, 1.0, time + dt,
                                      &AssemblerType::applyBCsBackwardEuler);
    std::function<VectorPtr(VectorPtr)> mFun =
          std::bind(&SteadySolver<AssemblerType>::assembleF,
                    this, time + dt, std::placeholders::_1, dt);

    std::function<GlobalBlockMatrix(VectorPtr)> mJac =
          std::bind(&SteadySolver<AssemblerType>::assembleJac,
                    this, time + dt, std::placeholders::_1, dt);

    std::function<GlobalBlockMatrix(VectorPtr)> mJacPrec =
          std::bind(&SteadySolver<AssemblerType>::assembleJacPrec,
                    this, time + dt, std::placeholders::_1, dt);

    double tol = M_datafile("newton_method/tol", 1e-5);
    double maxit = M_datafile("newton_method/maxit", 5);
    this->solveNonLinearSystem(mFun, mJac, M_solution, tol, maxit, &mJacPrec);

    printlog(MAGENTA, "done\n", M_verbose);
}

template <class AssemblerType>
typename SteadySolver<AssemblerType>::VectorPtr
SteadySolver<AssemblerType>::
assembleF(const double& time, VectorPtr tentativeSol, const double& dt)
{
    std::string msg("assembling F\n");
    printlog(GREEN, msg, M_verbose);
    M_globalAssembler->setTimeAndPrevSolution(time, tentativeSol);
    VectorPtr retF(new Vector(M_solution->map()));
    retF->zero();

    *retF += *(M_globalAssembler->computeF());
    M_globalAssembler->applyBCsVector(retF, 0.0, time,
                                      &AssemblerType::applyBCsBackwardEuler);
    return retF;
}

template <class AssemblerType>
GlobalBlockMatrix
SteadySolver<AssemblerType>::
assembleJac(const double& time, VectorPtr tentativeSol, const double& dt)
{
    std::string msg("assembling Jacobian\n");
    printlog(GREEN, msg, M_verbose);

    M_globalAssembler->setTimeAndPrevSolution(time, tentativeSol, false);
    double diagonalCoefficient = 1.0;
    GlobalBlockMatrix retJac = M_globalAssembler->getJacobianF(true,
                                                          &diagonalCoefficient);

    return retJac;
}

template <class AssemblerType>
GlobalBlockMatrix
SteadySolver<AssemblerType>::
assembleJacPrec(const double& time, VectorPtr tentativeSol, const double& dt)
{
    std::string msg("assembling preconditioner Jacobian\n");
    printlog(GREEN, msg, M_verbose);
    // we don't assemble the blocks again because this function is called
    // after assembleF (and the stabilization blocks are already assembled there)
    M_globalAssembler->setTimeAndPrevSolution(time, tentativeSol, false);
    double diagonalCoefficient = 1.0;
    GlobalBlockMatrix retJac = M_globalAssembler->getJacobianFprec(true,
                                                          &diagonalCoefficient);
    return retJac;
}

}  // namespace RedMA
