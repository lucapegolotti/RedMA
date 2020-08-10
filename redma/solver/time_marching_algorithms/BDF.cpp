#include "BDF.hpp"

namespace RedMA
{

BDF::
BDF(const DataContainer& data) :
  aTimeMarchingAlgorithm(data),
  M_order(data("time_discretization/order",2))
{
    M_useExtrapolation = this->M_data("time_discretization/use_extrapolation", 0);

    // if we set this we save the evaluation of the residual at the end of resolution
    // (important for rb method)
    if (M_useExtrapolation)
        this->M_systemSolver.isLinearProblem();
}

BDF::
BDF(const DataContainer& data, SHP(FunProvider) funProvider) :
  aTimeMarchingAlgorithm(data, funProvider),
  M_order(data("time_discretization/order",2))
{
    setup(this->M_funProvider->getZeroVector());

    M_useExtrapolation = this->M_data("time_discretization/use_extrapolation", 0);

    if (M_useExtrapolation)
        this->M_systemSolver.isLinearProblem();
}

void
BDF::
setup(const BlockVector& zeroVector)
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
        M_coefficients[0] = -1.0;
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

BlockVector
BDF::
computeExtrapolatedSolution()
{
    BlockVector extrapolatedSolution;
    extrapolatedSolution.hardCopy(M_prevSolutions[0]);

    if (M_order == 1)
    {

    }
    else if (M_order == 2)
    {
        extrapolatedSolution *= 2;
        extrapolatedSolution -= M_prevSolutions[1];
    }
    else if (M_order == 3)
    {
        extrapolatedSolution *= 3;
        extrapolatedSolution -= M_prevSolutions[1] * 3;
        extrapolatedSolution += M_prevSolutions[2];
    }
    else
        throw new Exception("BDF scheme of requested order not implemented");

    return extrapolatedSolution;
}

BlockVector
BDF::
advance(const double& time, double& dt, int& status)
{
    typedef BlockVector               BV;
    typedef BlockMatrix               BM;

    // we set the initial guess equal to the last solution
    // keep in mind that this MUST be a hard copy
    BV initialGuess = computeExtrapolatedSolution();
    // initialGuess.hardCopy(M_prevSolutions[0]);

    this->M_funProvider->applyDirichletBCs(time+dt, initialGuess);
    FunctionFunctor<BV,BV> fct(
        [this,time,dt](BV sol)
    {
        BM mass(this->M_funProvider->getMass(time+dt, sol));
        if (M_useExtrapolation)
            this->M_funProvider->setExtrapolatedSolution(computeExtrapolatedSolution());
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
        // the previous solution satisfies the boundary conditions so we search
        // for an increment with 0bcs
        this->M_funProvider->apply0DirichletBCs(retVec);

        return retVec;
    });

    FunctionFunctor<BV,BM> jac(
        [this,time,dt](BV sol)
    {
        // here the choice of hard copies is compulsory
        BM retMat;

        if (M_useExtrapolation)
            this->M_funProvider->setExtrapolatedSolution(computeExtrapolatedSolution());

        retMat.hardCopy(this->M_funProvider->getJacobianRightHandSide(time+dt, sol));
        retMat *= (-1. * M_rhsCoeff * dt);
        retMat += this->M_funProvider->getMass(time+dt, sol);
        retMat += this->M_funProvider->getMassJacobian(time+dt, sol);

        // this is the part relative to the previous steps
        unsigned int count = 0;
        for (BV vec : M_prevSolutions)
        {
            retMat += this->M_funProvider->getMassJacobian(time+dt, vec) * M_coefficients[count];
            count++;
        }

        return retMat;
    });

    BV sol = this->M_systemSolver.solve(fct, jac, initialGuess, status);

    if (status != 0)
        throw new Exception("Solver has not converged!");

    std::vector<SolverStatistics> statistics = this->M_systemSolver.getSolverStatistics();
    this->dumpSolverStatistics(statistics, time+dt);

    return sol;
}

BlockVector
BDF::
computeDerivative(const BlockVector& solnp1, double& dt)
{
    typedef BlockVector               BV;

    BlockVector retVec;
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

void
BDF::
shiftSolutions(const BlockVector& sol)
{
    // shift solutions
    std::vector<BlockVector> newPrevSolutions(M_order);
    newPrevSolutions[0].hardCopy(sol);

    for (unsigned int i = 0; i < M_order-1; i++)
        newPrevSolutions[i+1].softCopy(M_prevSolutions[i]);

    M_prevSolutions = newPrevSolutions;
}

}