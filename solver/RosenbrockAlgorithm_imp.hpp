// implementation of template class

namespace RedMA
{

template <class AssemblerType>
RosenbrockAlgorithm<AssemblerType>::
RosenbrockAlgorithm(const GetPot& datafile) :
  TimeMarchingAlgorithm<AssemblerType>(datafile),
  M_coefficients(datafile("time_discretization/scheme", "ROS3Pw"))
{
}

template <class AssemblerType>
void
RosenbrockAlgorithm<AssemblerType>::
solveTimestep(const double &time, double &dt,
              GlobalAssemblerType& assembler,
              const LinearSolver& linearSolver)
{
    typedef LifeV::VectorEpetra         VectorEpetra;
    unsigned int s = M_coefficients.numberStages();
    assembler.setTimeAndPrevSolution(time, M_solution);

    MapEpetraPtr globalMap = assembler.getGlobalMap();
    MatrixPtr globalMass = assembler.getGlobalMass();
    MatrixPtr globalJac = assembler.getJacobianF();

    MatrixPtr systemMatrix(new Matrix(*globalJac));
    *systemMatrix *= (-dt * M_coefficients.gamma());
    *systemMatrix += (*globalMass);

    std::vector<VectorPtr> stages(s);

    VectorPtr Fder = assembler.computeFder();

    for (int i = 0; i < s; i++)
    {
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
        // VectorPtr F = assembler.computeF(time + dt * alphai, yTilde);
        VectorPtr F = assembler.computeF();
        *F *= (M_coefficients.gamma() * dt);

        VectorPtr sumStages(new Vector(*globalMap));
        sumStages->zero();
        for (int j = 0; j < i; j++)
        {
            VectorEpetra part = M_coefficients.invGamma(i,j) * (*(stages[j]));
            *sumStages += (-M_coefficients.gamma() * part);
        }

        VectorEpetra prod = (*globalMass) * (*sumStages);
        *F += prod;

        double coeff = M_coefficients.gamma() * gammai * dt * dt;
        *Fder *= (coeff);

        *F += *Fder;

        *Fder *= (1.0/coeff) * (*Fder);

        // here we need to apply the bcs

        VectorPtr newStage(new Vector(*globalMap));
        // linearSolver.solve(newStage, systemMatrix, F);
        stages[i] = newStage;
    }

    // retrieve solution
    for (int i = 0; i < s; i++)
    {
        VectorEpetra part = M_coefficients.mHigh(i) * (*(stages[i]));
        *M_solution += (dt * part);
    }

}

}  // namespace RedMA
