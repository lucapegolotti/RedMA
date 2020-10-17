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
setup(const SHP(aVector)& zeroVector)
{
    M_coefficients.reserve(M_order);

    for (unsigned int i = 0; i < M_order; i++)
    {
        SHP(BlockVector) zeroVectorCopy(new BlockVector(0));
        zeroVectorCopy->hardCopy(zeroVector);
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

SHP(aVector)
BDF::
computeExtrapolatedSolution()
{
    SHP(BlockVector) extrapolatedSolution(new BlockVector(0));
    extrapolatedSolution->hardCopy(M_prevSolutions[0]);
    // std::cout << "extrapolatedSolution open?" << extrapolatedSolution->isOpen() << std::endl << std::flush;
    // std::cout << "extrapolatedSolution open?" << extrapolatedSolution->isOpen() << std::endl << std::flush;
    if (M_order == 1)
    {

    }
    else if (M_order == 2)
    {
        // std::cout << "extrapolatedSolution open?" << extrapolatedSolution->isOpen() << std::endl << std::flush;
        extrapolatedSolution->multiplyByScalar(2);
        M_prevSolutions[1]->multiplyByScalar(-1);
        extrapolatedSolution->add(M_prevSolutions[1]);
        M_prevSolutions[1]->multiplyByScalar(-1);
    }
    else if (M_order == 3)
    {
        extrapolatedSolution->multiplyByScalar(3);
        // M_prevSolutions[1]->close();
        M_prevSolutions[1]->multiplyByScalar(-3);
        extrapolatedSolution->add(M_prevSolutions[1]);
        M_prevSolutions[1]->multiplyByScalar(-1.0/3.0);
        // M_prevSolutions[2]->close();
        extrapolatedSolution->add(M_prevSolutions[2]);
    }
    else
        throw new Exception("BDF scheme of requested order not implemented");

    return extrapolatedSolution;
}

SHP(aVector)
BDF::
advance(const double& time, double& dt, int& status)
{
    typedef SHP(aVector)               BV;
    typedef SHP(aMatrix)               BM;

    // we set the initial guess equal to the last solution
    // keep in mind that this MUST be a hard copy
    BV initialGuess = computeExtrapolatedSolution();

    this->M_funProvider->applyDirichletBCs(time+dt, initialGuess);

    FunctionFunctor<BV,BV> fct(
        [this,time,dt](BV sol)
    {
        BM mass(this->M_funProvider->getMass(time+dt, sol));
        if (M_useExtrapolation)
            this->M_funProvider->setExtrapolatedSolution(computeExtrapolatedSolution());

        BV f(this->M_funProvider->getRightHandSide(time+dt, sol));

        BV prevContribution(new BlockVector(0));
        unsigned int count = 0;
        for (BV vec : M_prevSolutions)
        {
            BV vecCopy(new BlockVector(0));
            vecCopy->hardCopy(M_prevSolutions[count]);
            vecCopy->multiplyByScalar(M_coefficients[count]);
            prevContribution->add(vecCopy);
            count++;
        }
        prevContribution->add(sol);

        BV retVec(new BlockVector(0));
        retVec->hardCopy(mass->multiplyByVector(prevContribution));

        f->multiplyByScalar(-1. * M_rhsCoeff * dt);

        retVec->add(f);
        // the previous solution satisfies the boundary conditions so we search
        // for an increment with 0bcs
        this->M_funProvider->apply0DirichletBCs(retVec);

        return retVec;
    });

    FunctionFunctor<BV,BM> jac(
        [this,time,dt](BV sol)
    {
        // here the choice of hard copies is compulsory
        BM retMat(new BlockMatrix(1,1));

        if (M_useExtrapolation)
            this->M_funProvider->setExtrapolatedSolution(computeExtrapolatedSolution());

        retMat->hardCopy(this->M_funProvider->getJacobianRightHandSide(time+dt, sol));
        retMat->multiplyByScalar(-1. * M_rhsCoeff * dt);
        retMat->add(this->M_funProvider->getMass(time+dt, sol));
        retMat->add(this->M_funProvider->getMassJacobian(time+dt, sol));

        // this is the part relative to the previous steps
        unsigned int count = 0;
        for (BV vec : M_prevSolutions)
        {
            auto massjac = this->M_funProvider->getMassJacobian(time+dt, vec);
            massjac->multiplyByScalar(M_coefficients[count]);
            retMat->add(massjac);
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

SHP(aVector)
BDF::
computeDerivative(const SHP(aVector)& solnp1, double& dt)
{
    typedef SHP(BlockVector)               BV;

    BV retVec(new BlockVector(0));
    retVec->hardCopy(solnp1);

    unsigned int count = 0;
    for (BV vec : M_prevSolutions)
    {
        BV vecCopy(new BlockVector(0));
        vecCopy->hardCopy(vec);
        vecCopy->multiplyByScalar(M_coefficients[count]);
        retVec->add(vecCopy);
        count++;
    }

    retVec->multiplyByScalar(1.0/(dt * M_rhsCoeff));
    return retVec;
}

void
BDF::
shiftSolutions(const SHP(aVector)& sol)
{
    // shift solutions
    std::vector<SHP(BlockVector)> newPrevSolutions(M_order);
    newPrevSolutions[0].reset(static_cast<BlockVector*>(sol->clone()));

    for (unsigned int i = 0; i < M_order-1; i++)
        newPrevSolutions[i+1].reset(new BlockVector(*M_prevSolutions[i]));

    M_prevSolutions = newPrevSolutions;
}

}
