// implementation of template class

namespace RedMA
{

template <class AssemblerType>
BackwardEuler<AssemblerType>::
BackwardEuler(const GetPot& datafile,
              GlobalAssemblerType* assembler,
              commPtr_Type comm,
              bool verbose) :
  TimeMarchingAlgorithm<AssemblerType>(datafile, assembler, comm, verbose),
  M_coefficients(datafile("time_discretization/scheme", "ROS2"))
{
    double diagonalCoefficient = 1.0;
    // first we store the mass with no boundary conditions
    M_massMatrixNoBCs = assembler->assembleGlobalMass(false);

    assembler->assembleGlobalMass(false, &diagonalCoefficient);
    M_prevSolution.reset(new Vector(assembler->getGlobalMap()));
    M_prevSolution->zero();
}

template <class AssemblerType>
void
BackwardEuler<AssemblerType>::
solveTimestep(const double &time, double &dt)
{
    typedef LifeV::VectorEpetra         VectorEpetra;

    std::string msg("[BackwardEuler] solving, time = ");
    msg += std::to_string(time) + " ...\n";
    printlog(MAGENTA, msg, M_verbose);
    *M_prevSolution = *M_solution;
    M_globalAssembler->applyBCsVector(M_solution, 1.0, time + dt,
                                      &AssemblerType::applyBCsBackwardEuler);
    std::function<VectorPtr(VectorPtr)> mFun =
          std::bind(&BackwardEuler<AssemblerType>::assembleF,
                    this, time + dt, std::placeholders::_1, dt);

    std::function<MatrixPtr(VectorPtr)> mJac =
          std::bind(&BackwardEuler<AssemblerType>::assembleJac,
                    this, time + dt, std::placeholders::_1, dt);

    double tol = M_datafile("newton_method/tol", 1e-5);
    double maxit = M_datafile("newton_method/maxit", 5);
    this->solveNonLinearSystem(mFun, mJac, M_solution, tol, maxit);

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

    *retF = (*M_massMatrixNoBCs) * (*M_prevSolution);
    *retF *= (-1.0);
    *retF += (*M_massMatrixNoBCs) * (*tentativeSol);
    *retF -= (dt) * (*(M_globalAssembler->computeF()));
    M_globalAssembler->applyBCsVector(retF, 0.0, time,
                                      &AssemblerType::applyBCsBackwardEuler);
    return retF;
}

template <class AssemblerType>
typename BackwardEuler<AssemblerType>::MatrixPtr
BackwardEuler<AssemblerType>::
assembleJac(const double& time, VectorPtr tentativeSol, const double& dt)
{
    std::string msg("assembling Jacobian\n");
    printlog(GREEN, msg, M_verbose);
    M_globalAssembler->setTimeAndPrevSolution(time, tentativeSol);
    double diagonalCoefficient = 0;
    MatrixPtr retJac = M_globalAssembler->getJacobianF(true,
                                                       &diagonalCoefficient);
    *retJac *= (-dt);
    *retJac += *M_globalAssembler->getGlobalMass();

    retJac->spy("sytemMatrix");
    return retJac;
}

}  // namespace RedMA