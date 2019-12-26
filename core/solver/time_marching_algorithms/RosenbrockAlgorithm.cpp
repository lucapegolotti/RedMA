#include <RosenbrockAlgorithm.hpp>

namespace RedMA
{

RosenbrockAlgorithm::
RosenbrockAlgorithm(const GetPot& datafile,
                    GlobalAssemblerType* assembler,
                    commPtr_Type comm,
                    bool verbose) :
  TimeMarchingAlgorithm(datafile, assembler, comm, verbose),
  M_coefficients(datafile("time_discretization/scheme", "ROS2"))
{
    double diagonalCoefficient = 1.0;
    // first we store the mass with no boundary conditions
    M_massMatrixNoBCs = assembler->assembleGlobalMass(false);

    assembler->assembleGlobalMass(false, &diagonalCoefficient);

    this->M_order = M_coefficients.order();
}

void
RosenbrockAlgorithm::
solveTimestep(const double &time, double &dt)
{
    typedef LifeV::VectorEpetra         VectorEpetra;

    std::string msg("[RosenbrockAlgorithm] solving, time = ");
    msg += std::to_string(time) + " ...\n";
    printlog(MAGENTA, msg, M_verbose);
    unsigned int s = M_coefficients.numberStages();
    M_globalAssembler->setTimeAndPrevSolution(time, M_solution);

    MapEpetraPtr globalMap = M_globalAssembler->getGlobalMap();
    BlockMatrix globalMass = M_globalAssembler->getGlobalMass();

    double diagonalCoefficient = 0.0;
    BlockMatrix globalJac = M_globalAssembler->getJacobianF(true,
                                                          &diagonalCoefficient);
    BlockMatrix systemMatrix(globalJac);
    systemMatrix *= (-dt * M_coefficients.gamma());
    // systemMatrix->openCrsMatrix();
    systemMatrix.add(globalMass);
    // systemMatrix->globalAssemble();
    std::vector<VectorPtr> stages(s);

    VectorPtr Fder = M_globalAssembler->computeFder_();

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
        VectorPtr F = M_globalAssembler->computeF_();
        *F *= (M_coefficients.gamma() * dt);

        VectorPtr sumStages(new Vector(*globalMap));
        sumStages->zero();
        for (int j = 0; j < i; j++)
        {
            VectorEpetra part = M_coefficients.invGamma(i,j) * (*(stages[j]));
            *sumStages += (-M_coefficients.gamma() * part);
        }

        // VectorEpetra prod = (M_massMatrixNoBCs) * (*sumStages);
        VectorEpetra prod = M_globalAssembler->getGlobalMass() * (*sumStages);
        *F += prod;

        double coeff = M_coefficients.gamma() * gammai * dt * dt;
        *Fder *= coeff;

        *F += *Fder;

        *Fder *= (1.0/coeff); // * (*Fder);

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
