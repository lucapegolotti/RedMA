namespace RedMA
{

template <class InVectorType, class InMatrixType>
BDF<InVectorType, InMatrixType>::
BDF(const DataContainer& data) :
  aTimeMarchingAlgorithm<InVectorType, InMatrixType>(data),
  M_order(data("time_discretization/order",2))
{
}

template <class InVectorType, class InMatrixType>
BDF<InVectorType, InMatrixType>::
BDF(const DataContainer& data, SHP(FunProvider) funProvider) :
  aTimeMarchingAlgorithm<InVectorType, InMatrixType>(data, funProvider),
  M_order(data("time_discretization/order",2))
{
    setup(this->M_funProvider->getZeroVector());
}

template <class InVectorType, class InMatrixType>
void
BDF<InVectorType, InMatrixType>::
setup(const BlockVector<InVectorType>& zeroVector)
{
    M_coefficients.reserve(M_order);

    for (unsigned int i = 0; i < M_order; i++)
    {
        BlockVector<InVectorType> zeroVectorCopy;
        zeroVectorCopy.hardCopy(zeroVector);
        M_prevSolutions.push_back(zeroVectorCopy);
    }

    if (M_order == 1)
    {
        M_coefficients[0] = 1.0;
        M_rhsCoeff = 1.0;
    }
    else if (M_order == 2)
    {
        M_coefficients[0] = -4.0/3.0;
        M_coefficients[1] = 1.0/3.0;
        M_rhsCoeff = 2.0/3.0;
    }
    else if (M_order == 3)
    {
        M_coefficients[0] = -18.0/11.0;
        M_coefficients[1] = 9.0/11.0;
        M_coefficients[2] = -2.0/11.0;
        M_rhsCoeff = 6.0/11.0;
    }
    else
    {
        throw new Exception("BDF scheme of requested order not implemented");
    }
}

template <class InVectorType, class InMatrixType>
BlockVector<InVectorType>
BDF<InVectorType, InMatrixType>::
advance(const double& time, double& dt, int& status)
{
    typedef BlockVector<InVectorType>               BV;
    typedef BlockMatrix<InMatrixType>               BM;

    // we set the initial guess equal to the last solution
    // keep in mind that this MUST be a hard copy
    BV initialGuess;
    initialGuess.hardCopy(M_prevSolutions[0]);

    FunctionFunctor<BV,BV> fct(
        [this,time,dt](BV sol)
    {
        BM mass(this->M_funProvider->getMass(time+dt, sol));
        BV f(this->M_funProvider->getRightHandSide(time+dt, sol));
        BV prevContribution;

        unsigned int count = 0;
        for (BV vec : M_prevSolutions)
        {
            prevContribution += vec * M_coefficients[count];
            count++;
        }

        BV retVec;
        retVec.softCopy(mass * (sol + prevContribution));
        f *= (-1. * M_rhsCoeff * dt);
        retVec += f;

        return retVec;
    });

    FunctionFunctor<BV,BM> jac(
        [this,time,dt](BV sol)
    {
        // here the choice of hard copies is compulsory
        BM retMat;

        retMat.hardCopy(this->M_funProvider->getJacobianRightHandSide(time+dt, sol));
        retMat *= (-1. * M_rhsCoeff * dt);
        retMat += this->M_funProvider->getMass(time+dt, sol);
        retMat += this->M_funProvider->getMassJacobian(time+dt, sol);

        return retMat;
    });

    BV sol = this->M_systemSolver.solve(fct, jac, initialGuess, status);

    if (status != 0)
        throw new Exception("Solver has not converged!");

    std::vector<SolverStatistics> statistics = this->M_systemSolver.getSolverStatistics();
    this->dumpSolverStatistics(statistics, time+dt);

    return sol;
}

template <class InVectorType, class InMatrixType>
BlockVector<InVectorType>
BDF<InVectorType, InMatrixType>::
computeDerivative(const BlockVector<InVectorType>& solnp1, double& dt)
{
    typedef BlockVector<InVectorType>               BV;

    BlockVector<InVectorType> retVec;
    retVec.hardCopy(solnp1);

    unsigned int count = 0;
    for (BV vec : M_prevSolutions)
    {
        retVec += vec * M_coefficients[count];
        count++;
    }

    retVec *= (1.0/(dt * M_rhsCoeff));
    return retVec;
}

template <class InVectorType, class InMatrixType>
void
BDF<InVectorType, InMatrixType>::
shiftSolutions(const BlockVector<InVectorType>& sol)
{
    // shift solutions
    std::vector<BlockVector<InVectorType>> newPrevSolutions(M_order);
    newPrevSolutions[0].hardCopy(sol);

    for (unsigned int i = 0; i < M_order-1; i++)
        newPrevSolutions[i+1].softCopy(M_prevSolutions[i]);

    M_prevSolutions = newPrevSolutions;
}

}
