// implementation of template class

namespace RedMA
{

template <class AssemblerType>
RosenbrockAlgorithm<AssemblerType>::
RosenbrockAlgorithm(const GetPot& datafile,
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
}

template <class AssemblerType>
void
RosenbrockAlgorithm<AssemblerType>::
solveTimestep(const double &time, double &dt)
{
    typedef LifeV::VectorEpetra         VectorEpetra;

    std::string msg("[RosenbrockAlgorithm] solving, time = ");
    msg += std::to_string(time) + " ...\n";
    printlog(MAGENTA, msg, M_verbose);
    unsigned int s = M_coefficients.numberStages();
    M_globalAssembler->setTimeAndPrevSolution(time, M_solution);

    MapEpetraPtr globalMap = M_globalAssembler->getGlobalMap();
    MatrixPtr globalMass = M_globalAssembler->getGlobalMass();

    double diagonalCoefficient = 0.0;
    MatrixPtr globalJac = M_globalAssembler->getJacobianF(true, &diagonalCoefficient);
    MatrixPtr systemMatrix(new Matrix(*globalJac));
    *systemMatrix *= (-dt * M_coefficients.gamma());
    // systemMatrix->openCrsMatrix();
    *systemMatrix += (*globalMass);
    systemMatrix->spy("systemMatrix");
    // systemMatrix->globalAssemble();
    std::vector<VectorPtr> stages(s);

    VectorPtr Fder = M_globalAssembler->computeFder();

    for (int i = 0; i < s; i++)
    {
        std::string msg_("Starting stage #");
        msg_ += std::to_string(i+1) + "\n";
        printlog(GREEN, msg_, M_verbose);
        VectorPtr yTilde(new Vector(*globalMap));
        *yTilde = *M_solution;

        double alphai = 0;
        double gammai = 0;

        for (int j = 0; j <= i; j++)
        {
            alphai += M_coefficients.alpha(i,j);
            gammai += M_coefficients.gamma(i,j);
        }

        for (int j = 0; j < i; j++)
        {
            VectorEpetra part = M_coefficients.alphaHat(i,j) * (*(stages[j]));
            *yTilde += part;
        }

        M_globalAssembler->setTimeAndPrevSolution(time + dt * alphai, yTilde);
        VectorPtr F = M_globalAssembler->computeF();
        *F *= (M_coefficients.gamma() * dt);

        VectorPtr sumStages(new Vector(*globalMap));
        sumStages->zero();
        for (int j = 0; j < i; j++)
        {
            VectorEpetra part = M_coefficients.invGamma(i,j) * (*(stages[j]));
            *sumStages += (-M_coefficients.gamma() * part);
        }

        VectorEpetra prod = (*M_massMatrixNoBCs) * (*sumStages);
        *F += prod;

        double coeff = M_coefficients.gamma() * gammai * dt * dt;
        *Fder *= coeff;

        *F += *Fder;

        *Fder *= (1.0/coeff) * (*Fder);

        // here we need to apply the bcs to the right hand side
        M_globalAssembler->applyBCsRhsRosenbrock(F, yTilde, time, dt,
                                                 alphai, gammai);
        VectorPtr newStage(new Vector(*globalMap));
        solveLinearSystem(systemMatrix, F, newStage);
        stages[i] = newStage;
    }
    // retrieve solution
    for (int i = 0; i < s; i++)
    {
        VectorEpetra part = M_coefficients.mHigh(i) * (*(stages[i]));
        *M_solution += part;
    }
    printlog(MAGENTA, "done\n", M_verbose);
}

}  // namespace RedMA
