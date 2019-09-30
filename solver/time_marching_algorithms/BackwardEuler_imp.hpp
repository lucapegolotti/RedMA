// implementation of template class

namespace RedMA
{

template <class AssemblerType>
BackwardEuler<AssemblerType>::
BackwardEuler(const GetPot& datafile,
              GlobalAssemblerType* assembler,
              commPtr_Type comm,
              bool verbose) :
  TimeMarchingAlgorithm<AssemblerType>(datafile, assembler, comm, verbose)
{
    double diagonalCoefficient = 1.0;
    // first we store the mass with no boundary conditions
    assembler->assembleGlobalMass(false, &diagonalCoefficient);
    M_prevSolution.reset(new Vector(assembler->getGlobalMap()));
    M_prevSolution->zero();

    this->M_order = 1;
}

template <class AssemblerType>
void
BackwardEuler<AssemblerType>::
solveTimestep(const double &time, double &dt)
{
    std::string msg("[BackwardEuler] solving, time = ");
    msg += std::to_string(time) + " ...\n";
    printlog(MAGENTA, msg, M_verbose);
    *M_prevSolution = *M_solution;
    M_globalAssembler->applyBCsVector(M_solution, 1.0, time + dt,
                                      &AssemblerType::applyBCsBackwardEuler);
    std::function<VectorPtr(VectorPtr)> mFun =
          std::bind(&BackwardEuler<AssemblerType>::assembleF,
                    this, time + dt, std::placeholders::_1, dt);

    std::function<GlobalBlockMatrix(VectorPtr)> mJac =
          std::bind(&BackwardEuler<AssemblerType>::assembleJac,
                    this, time + dt, std::placeholders::_1, dt);

    double tol = M_datafile("newton_method/tol", 1e-5);
    double maxit = M_datafile("newton_method/maxit", 5);
    this->solveNonLinearSystem(mFun, mJac, M_solution, tol, maxit);

    // M_globalAssembler->checkResidual(M_solution, M_prevSolution, dt);

    printlog(MAGENTA, "done\n", M_verbose);
}

template <class AssemblerType>
typename BackwardEuler<AssemblerType>::VectorPtr
BackwardEuler<AssemblerType>::
assembleF(const double& time, VectorPtr tentativeSol, const double& dt)
{
    std::string msg("assembling F\n");
    printlog(GREEN, msg, M_verbose);
    M_globalAssembler->setTimeAndPrevSolution(time, tentativeSol);
    VectorPtr retF(new Vector(M_solution->map()));
    retF->zero();

    *retF = (M_globalAssembler->getGlobalMass()) * (*M_prevSolution);
    *retF *= (-1.0);
    *retF += (M_globalAssembler->getGlobalMass()) * (*tentativeSol);

    *retF -= (dt) * (*(M_globalAssembler->computeF()));
    M_globalAssembler->applyBCsVector(retF, 0.0, time,
                                      &AssemblerType::applyBCsBackwardEuler);

    return retF;
}

template <class AssemblerType>
GlobalBlockMatrix
BackwardEuler<AssemblerType>::
assembleJac(const double& time, VectorPtr tentativeSol, const double& dt)
{
    std::string msg("assembling Jacobian\n");
    printlog(GREEN, msg, M_verbose);
    // we don't assemble the blocks again because this function is called
    // after assembleF (and the stabilization blocks are already assembled there)
    M_globalAssembler->setTimeAndPrevSolution(time, tentativeSol, false);
    double diagonalCoefficient = 0;
    GlobalBlockMatrix retJac = M_globalAssembler->getJacobianF(true,
                                                          &diagonalCoefficient);

    retJac *= (-dt);

    // trick to compute the jacobian of the part corresponding to H(un+1) * un
    M_globalAssembler->setTimeAndPrevSolution(time, M_prevSolution, false);
    GlobalBlockMatrix massUn = M_globalAssembler->getGlobalMassJacVelocity();

    massUn *= (-1);

    retJac.add(massUn);
    retJac.add(M_globalAssembler->getGlobalMassJac());
    return retJac;
}

}  // namespace RedMA
