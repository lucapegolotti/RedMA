#include "GeneralizedAlphaMethod.hpp"

namespace RedMA
{

GeneralizedAlphaMethod::
GeneralizedAlphaMethod(const DataContainer& data) :
  aTimeMarchingAlgorithm(data),
  M_order(2)
{

}

GeneralizedAlphaMethod::
GeneralizedAlphaMethod(const DataContainer& data, SHP(FunProvider) funProvider) :
  aTimeMarchingAlgorithm(data, funProvider),
  M_order(2)
{
  setup(this->M_funProvider->getZeroVector());
}

void
GeneralizedAlphaMethod::
setup(const SHP(aVector)& zeroVector)
{
    // M_prevSolution.hardCopy(zeroVector);
    // M_prevDerivative.hardCopy(zeroVector);
    //
    // M_rhoinf = this->M_data("time_discretization/rhoinf", 0.5);
    //
    // M_alpham = 0.5 * (3.0 - M_rhoinf) / (1.0 + M_rhoinf);
    // M_alphaf = 1.0 / (1.0 + M_rhoinf);
    // M_gamma = 0.5 + M_alpham - M_alphaf;
}

SHP(aVector)
GeneralizedAlphaMethod::
advance(const double& time, double& dt, int& status)
{
    // typedef BlockVector               BV;
    // typedef BlockMatrix               BM;
    //
    // // we are actually solving for the derivative
    // BV initialGuess;
    // initialGuess.hardCopy(M_prevDerivative);
    //
    // std::cout << "think about boundary conditions for generalized alpha " << std::endl << std::flush;
    // exit(1);
    // FunctionFunctor<BV,BV> fct(
    //     [this,time,dt](BV dersol)
    // {
    //     BV solnp1 = computesolnp1(dersol, dt);
    //     BV solnpalphaf = computesolnpalphaf(solnp1);
    //     BV dersolnpalpham = computedersolnpalpham(dersol);
    //
    //     BM mass(this->M_funProvider->getMass(time + dt*M_alphaf, solnpalphaf));
    //
    //     // attention: the neumann term could actually be evaluated at time + dt*M_alphaf
    //     // ask ju about this
    //     BV f(this->M_funProvider->getRightHandSide(time + dt*M_alphaf, solnpalphaf));
    //
    //     BV retVec = mass * dersolnpalpham;
    //     retVec -= f;
    //
    //     return retVec;
    // });
    //
    // FunctionFunctor<BV,BM> jac(
    //     [this,time,dt](BV dersol)
    // {
    //     // here the choice of hard copies is compulsory
    //     BV solnp1 = computesolnp1(dersol, dt);
    //     BV solnpalphaf = computesolnpalphaf(solnp1);
    //     BV dersolnpalpham = computedersolnpalpham(dersol);
    //
    //     // here the choice of hard copies is compulsory
    //     BM retMat;
    //
    //     retMat.hardCopy(this->M_funProvider->getJacobianRightHandSide(time + dt*M_alphaf, solnpalphaf));
    //     // we are computing the jacobian actually with respect to dersol -> we multiply
    //     // by the derivative of solnp1 wrt dersol
    //     const double coeff = M_gamma * M_alphaf * dt;
    //     retMat *= (-1. * coeff);
    //     retMat += this->M_funProvider->getMass(time + dt*M_alphaf, solnpalphaf) * M_alpham;
    //     retMat += (this->M_funProvider->getMassJacobian(time + dt*M_alphaf, solnpalphaf) * coeff);
    //
    //     return retMat;
    // });
    //
    // BV dersol = this->M_systemSolver.solve(fct, jac, initialGuess, status);
    //
    // if (status != 0)
    //     throw new Exception("Solver has not converged!");
    //
    // BV sol = computesolnp1(dersol, dt);
    //
    // std::vector<SolverStatistics> statistics = this->M_systemSolver.getSolverStatistics();
    // this->dumpSolverStatistics(statistics, time+dt);
    //
    // return sol;
}

SHP(aVector)
GeneralizedAlphaMethod::
computesolnp1(SHP(aVector) dersol, const double& dt)
{
    // BlockVector solnp1;
    // solnp1.hardCopy(M_prevSolution);
    // solnp1 += (M_prevDerivative * dt);
    // solnp1 += (dersol - M_prevDerivative) * (M_gamma * dt);
    //
    // return solnp1;
}

SHP(aVector)
GeneralizedAlphaMethod::
computesolnpalphaf(SHP(aVector) solnp1)
{
    // BlockVector solnpalphaf;
    // solnpalphaf.hardCopy(M_prevSolution);
    // solnpalphaf += (solnp1 - M_prevSolution) * M_alphaf;
    //
    // return solnpalphaf;
}

SHP(aVector)
GeneralizedAlphaMethod::
computedersolnpalpham(SHP(aVector) dersol)
{
    // BlockVector dersolnpalpham;
    // dersolnpalpham.hardCopy(M_prevDerivative);
    // dersolnpalpham += (dersol - M_prevDerivative) * M_alpham;
    //
    // return dersolnpalpham;
}

SHP(aVector)
GeneralizedAlphaMethod::
computeDerivative(const SHP(aVector)& solnp1, double& dt)
{
    // BlockVector dersolnpalpham;
    // dersolnpalpham.hardCopy(M_prevDerivative);
    // dersolnpalpham += (solnp1 - M_prevSolution - (M_prevDerivative * dt)) * (1.0 / (M_gamma * dt));
    // return dersolnpalpham;
}

void
GeneralizedAlphaMethod::
shiftSolutions(const SHP(aVector)& sol)
{
    // double dt = this->M_data("time_discretization/dt", 0.01);
    // M_prevDerivative = computeDerivative(sol, dt);
    // M_prevSolution = sol;
}

}
