#include "BDF.hpp"

namespace RedMA
{

BDF::
BDF(const DataContainer& data) :
  aTimeMarchingAlgorithm(data),
  M_order(data("time_discretization/order",2)),
  M_extrapolationOrder(data("time_discretization/extrapolation_order", 2))
{
    M_useExtrapolation = this->M_data("time_discretization/use_extrapolation", 0);
    // if we set this we save the evaluation of the residual at the end of resolution
    // (important for rb method)
    if (M_useExtrapolation)
        this->setLinearSolver();
}

BDF::
BDF(const DataContainer& data, shp<FunProvider> funProvider) :
  aTimeMarchingAlgorithm(data, funProvider),
  M_order(data("time_discretization/order",2)),
  M_extrapolationOrder(data("time_discretization/extrapolation_order", 2))
{
    setup(this->M_funProvider->getZeroVector());

    M_useExtrapolation = this->M_data("time_discretization/use_extrapolation", 0);

    if (M_useExtrapolation)
        this->M_systemSolver.isLinearProblem();
}

BDF::
BDF(const DataContainer& data, const shp<aVector>& zeroVector):
  aTimeMarchingAlgorithm(data, zeroVector),
  M_order(data("time_discretization/order",2)),
  M_extrapolationOrder(data("time_discretization/extrapolation_order", 2))
{
    setup(zeroVector);

    M_useExtrapolation = this->M_data("time_discretization/use_extrapolation", 0);
    // if we set this we save the evaluation of the residual at the end of resolution
    // (important for rb method)
    if (M_useExtrapolation)
        this->M_systemSolver.isLinearProblem();
}

void
BDF::
setup(const shp<aVector>& zeroVector)
{
    for (unsigned int i = 0; i < std::max(M_order, M_extrapolationOrder); i++)
    {
        shp<BlockVector> zeroVectorCopy(new BlockVector(0));
        zeroVectorCopy->deepCopy(zeroVector);
        M_prevSolutions.push_back(zeroVectorCopy);
    }

    this->setBDFCoefficients();
    this->setExtrapolationCoefficients();

}

void
BDF::
setBDFCoefficients()
{
    M_coefficients.reserve(M_order);

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

void
BDF::
setExtrapolationCoefficients()
{
    M_extrapolationCoefficients.reserve(M_extrapolationOrder);

    if (M_extrapolationOrder == 1)
    {
        M_extrapolationCoefficients[0] = 1.0;
    }
    else if (M_extrapolationOrder == 2)
    {
        M_extrapolationCoefficients[0] = 2.0;
        M_extrapolationCoefficients[1] = -1.0;
    }
    else if (M_extrapolationOrder == 3)
    {
        M_extrapolationCoefficients[0] = 3.0;
        M_extrapolationCoefficients[1] = -3.0;
        M_extrapolationCoefficients[2] = 1.0;
    }
    else
    {
        throw new Exception("Extrapolation scheme of requested order not implemented");
    }
}

shp<aVector>
BDF::
computeExtrapolatedSolution() {

    if (M_extrapolationOrder <= 0 || M_extrapolationOrder > 3)
        throw new Exception("Extrapolation scheme of requested order not implemented");

    shp<BlockVector> extrapolatedSolution(new BlockVector(0));
    extrapolatedSolution->deepCopy(M_prevSolutions[0]);
    extrapolatedSolution->multiplyByScalar(M_extrapolationCoefficients[0]);

    for (unsigned int i = 1; i < M_extrapolationOrder; i++) {
        M_prevSolutions[i]->multiplyByScalar(M_extrapolationCoefficients[i]);
        extrapolatedSolution->add(M_prevSolutions[i]);
        M_prevSolutions[i]->multiplyByScalar(1.0 / M_extrapolationCoefficients[i]);
    }

    return extrapolatedSolution;

}

shp<aVector>
BDF::
combineOldSolutions()
{
    shp<BlockVector> oldSteps(new BlockVector(0));

    oldSteps->deepCopy(M_prevSolutions[0]);
    oldSteps->multiplyByScalar(M_coefficients[0]);

    for (unsigned int i = 1; i < M_order; i++) {
        M_prevSolutions[i]->multiplyByScalar(M_coefficients[i]);
        oldSteps->add(M_prevSolutions[i]);
        M_prevSolutions[i]->multiplyByScalar(1.0 / M_coefficients[i]);
    }

    return oldSteps;
}

shp<aVector>
BDF::
advance(const double& time, double& dt, int& status)
{
    typedef shp<aVector>               BV;
    typedef shp<aMatrix>               BM;

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
        for (unsigned int count = 0; count < M_order; count++)
        {
            BV vecCopy(new BlockVector(0));
            vecCopy->deepCopy(M_prevSolutions[count]);
            vecCopy->multiplyByScalar(M_coefficients[count]);
            prevContribution->add(vecCopy);
        }
        prevContribution->add(sol);

        BV retVec(new BlockVector(0));
        retVec->deepCopy(mass->multiplyByVector(prevContribution));

        f->multiplyByScalar(-1. * M_rhsCoeff * dt);
        retVec->add(f);
        // the previous solution satisfies the boundary conditions, so we search
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

        retMat->deepCopy(this->M_funProvider->getJacobianRightHandSide(time+dt, sol));
        retMat->multiplyByScalar(-1. * M_rhsCoeff * dt);
        retMat->add(this->M_funProvider->getMass(time+dt, sol));
        retMat->add(this->M_funProvider->getMassJacobian(time+dt, sol));

        // this is the part relative to the previous steps
        for (unsigned int count = 0; count < M_order; count++)
        {
            auto massjac = this->M_funProvider->getMassJacobian(time+dt, M_prevSolutions[count]);
            massjac->multiplyByScalar(M_coefficients[count]);
            retMat->add(massjac);
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

shp<aVector>
BDF::
simpleAdvance(const double &dt, const shp<BlockVector> &sol)
{
    if (M_order <= 0 || M_order > 3)
        throw new Exception("BDF scheme of requested order not implemented");

    shp<BlockVector> retVec(new BlockVector(2));

    // I update only field 0, as for displacements no pressure is defined (i.e. it equals 0)
    retVec->deepCopy(sol);
    retVec->multiplyByScalar(dt * M_rhsCoeff);

    shp<BlockVector> oldSteps(new BlockVector(2));
    oldSteps->deepCopy(M_prevSolutions[0]);
    oldSteps->block(0)->multiplyByScalar(M_coefficients[0]);

    for (unsigned int i = 1; i < M_order; i++) {
        M_prevSolutions[i]->multiplyByScalar(M_coefficients[i]);
        oldSteps->block(0)->add(M_prevSolutions[i]->block(0));
        M_prevSolutions[i]->multiplyByScalar(1.0 / M_coefficients[i]);
    }
    oldSteps->block(0)->multiplyByScalar(-1.0);

    retVec->block(0)->add(oldSteps->block(0));

    return retVec;
}

shp<aVector>
BDF::
computeDerivative(const shp<aVector>& solnp1, double& dt)
{
    typedef shp<BlockVector>               BV;

    BV retVec(new BlockVector(0));
    retVec->deepCopy(solnp1);

    for (unsigned int count = 0; count < M_order; count++)
    {
        BV vecCopy(new BlockVector(0));
        vecCopy->deepCopy(M_prevSolutions[count]);
        vecCopy->multiplyByScalar(M_coefficients[count]);
        retVec->add(vecCopy);
        count++;
    }

    retVec->multiplyByScalar(1.0/(dt * M_rhsCoeff));
    return retVec;
}

void
BDF::
shiftSolutions(const shp<aVector>& sol)
{
    // shift solutions
    std::vector<shp<BlockVector>> newPrevSolutions(std::max(M_order, M_extrapolationOrder));
    shp<BlockVector> newSol = dpcast<BlockVector>(sol);
    newPrevSolutions[0].reset(newSol->clone());

    for (unsigned int i = 0; i < std::max(M_order, M_extrapolationOrder)-1; i++)
        newPrevSolutions[i+1].reset(new BlockVector(*M_prevSolutions[i]));

    M_prevSolutions = newPrevSolutions;
}

std::vector<double>
BDF::
getCoefficients() const
{
    std::vector<double> retVec = M_coefficients;
    retVec.push_back(M_rhsCoeff);

    return retVec;
}

}
