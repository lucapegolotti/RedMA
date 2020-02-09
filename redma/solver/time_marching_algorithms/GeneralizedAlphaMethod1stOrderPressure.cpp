#include "GeneralizedAlphaMethod1stOrderPressure.hpp"

namespace RedMA
{

template <>
BlockVector<BlockVector<VectorEp>>
GeneralizedAlphaMethod1stOrderPressure<BlockVector<VectorEp>, BlockMatrix<MatrixEp>>::
computesolnpalphaf(BlockVector<BlockVector<VectorEp>> solnp1)
{
    BlockVector<BlockVector<VectorEp>> solnpalphaf;
    solnpalphaf.hardCopy(M_prevSolution);
    solnpalphaf += (solnp1 - M_prevSolution) * M_alphaf;

    for (unsigned int i = 0; i < M_nPrimalBlocks; i++)
        solnpalphaf.block(i).block(1).hardCopy(solnp1.block(i).block(1));

    return solnpalphaf;
}

template <>
BlockVector<BlockVector<VectorEp>>
GeneralizedAlphaMethod1stOrderPressure<BlockVector<VectorEp>, BlockMatrix<MatrixEp>>::
advance(const double& time, double& dt, int& status)
{
    typedef BlockVector<BlockVector<VectorEp>>               BV;
    typedef BlockMatrix<BlockMatrix<MatrixEp>>               BM;

    // we are actually solving for the derivative
    BV initialGuess;
    initialGuess.hardCopy(M_prevDerivative);

    FunctionFunctor<BV,BV> fct(
        [this,time,dt](BV dersol)
    {
        BM mass(this->M_funProvider->getMass(time + dt*M_alphaf, M_prevSolution));
        // we count number of primal blocks

        M_nPrimalBlocks = 0;
        while (M_nPrimalBlocks < mass.nRows() && mass.block(M_nPrimalBlocks,M_nPrimalBlocks).nRows() > 0)
            M_nPrimalBlocks++;

        BV solnp1 = computesolnp1(dersol, dt);
        BV solnpalphaf = computesolnpalphaf(solnp1);
        BV dersolnpalpham = computedersolnpalpham(dersol);

        mass.softCopy(this->M_funProvider->getMass(time + dt*M_alphaf, solnpalphaf));
        // attention: the neumann term could actually be evaluated at time + dt*M_alphaf
        // ask ju about this
        BV f(this->M_funProvider->getRightHandSide(time + dt*M_alphaf, solnpalphaf));

        BV retVec = mass * dersolnpalpham;
        retVec -= f;

        return retVec;
    });

    FunctionFunctor<BV,BM> jac(
        [this,time,dt](BV dersol)
    {
        // here the choice of hard copies is compulsory
        BV solnp1 = computesolnp1(dersol, dt);
        BV solnpalphaf = computesolnpalphaf(solnp1);
        BV dersolnpalpham = computedersolnpalpham(dersol);

        // here the choice of hard copies is compulsory
        BM retMat;

        retMat.hardCopy(this->M_funProvider->getJacobianRightHandSide(time + dt*M_alphaf, solnpalphaf));
        // we are computing the jacobian actually with respect to dersol -> we multiply
        // by the derivative of solnp1 wrt dersol
        const double coeff = M_gamma * M_alphaf * dt;
        retMat *= (-1. * coeff);

        retMat += (this->M_funProvider->getMassJacobian(time + dt*M_alphaf, solnpalphaf) * coeff);

        for (unsigned int i = 0; i < M_nPrimalBlocks; i++)
        {
            retMat.block(i,i).block(0,1) *= (1.0/M_alphaf);
            retMat.block(i,i).block(1,1) *= (1.0/M_alphaf);
        }

        retMat += this->M_funProvider->getMass(time + dt*M_alphaf, solnpalphaf) * M_alpham;

        return retMat;
    });

    BV dersol = this->M_systemSolver.solve(fct, jac, initialGuess, status);

    if (status != 0)
        throw new Exception("Solver has not converged!");

    BV sol = computesolnp1(dersol, dt);

    std::vector<SolverStatistics> statistics = this->M_systemSolver.getSolverStatistics();
    this->dumpSolverStatistics(statistics, time+dt);

    return sol;
}

template <>
BlockVector<BlockVector<VectorEp>>
GeneralizedAlphaMethod1stOrderPressure<BlockVector<VectorEp>, BlockMatrix<MatrixEp>>::
computeDerivative(const BlockVector<BlockVector<VectorEp>>& solnp1, double& dt)
{
    BlockVector<BlockVector<VectorEp>> dersolnpalpham;
    dersolnpalpham.hardCopy(M_prevDerivative);
    dersolnpalpham += (solnp1 - M_prevSolution - (M_prevDerivative * dt)) * (1.0 / (M_gamma * dt));

    for (unsigned int i = 0; i < M_nPrimalBlocks; i++)
    {
        dersolnpalpham.block(i).block(1) = (solnp1.block(i).block(1) - M_prevSolution.block(i).block(1));
        dersolnpalpham.block(i).block(1) *= (1.0/dt);
    }

    return dersolnpalpham;
}

template <>
BlockVector<VectorEp>
GeneralizedAlphaMethod1stOrderPressure<VectorEp, MatrixEp>::
computesolnpalphaf(BlockVector<VectorEp> solnp1)
{
    BlockVector<VectorEp> solnpalphaf;
    solnpalphaf.hardCopy(M_prevSolution);
    solnpalphaf += (solnp1 - M_prevSolution) * M_alphaf;

    solnpalphaf.block(1).hardCopy(solnp1.block(1));

    return solnpalphaf;
}

template <>
BlockVector<VectorEp>
GeneralizedAlphaMethod1stOrderPressure<VectorEp, MatrixEp>::
advance(const double& time, double& dt, int& status)
{
    typedef BlockVector<VectorEp>               BV;
    typedef BlockMatrix<MatrixEp>               BM;

    // we are actually solving for the derivative
    BV initialGuess;
    initialGuess.hardCopy(M_prevDerivative);

    FunctionFunctor<BV,BV> fct(
        [this,time,dt](BV dersol)
    {
        BV solnp1 = computesolnp1(dersol, dt);
        BV solnpalphaf = computesolnpalphaf(solnp1);
        BV dersolnpalpham = computedersolnpalpham(dersol);

        BM mass(this->M_funProvider->getMass(time + dt*M_alphaf, solnpalphaf));
        // attention: the neumann term could actually be evaluated at time + dt*M_alphaf
        // ask ju about this
        BV f(this->M_funProvider->getRightHandSide(time + dt*M_alphaf, solnpalphaf));

        BV retVec = mass * dersolnpalpham;
        retVec -= f;

        return retVec;
    });

    FunctionFunctor<BV,BM> jac(
        [this,time,dt](BV dersol)
    {
        // here the choice of hard copies is compulsory
        BV solnp1 = computesolnp1(dersol, dt);
        BV solnpalphaf = computesolnpalphaf(solnp1);
        BV dersolnpalpham = computedersolnpalpham(dersol);

        // here the choice of hard copies is compulsory
        BM retMat;

        retMat.hardCopy(this->M_funProvider->getJacobianRightHandSide(time + dt*M_alphaf, solnpalphaf));
        // we are computing the jacobian actually with respect to dersol -> we multiply
        // by the derivative of solnp1 wrt dersol
        const double coeff = M_gamma * M_alphaf * dt;
        retMat *= (-1. * coeff);

        retMat += (this->M_funProvider->getMassJacobian(time + dt*M_alphaf, solnpalphaf) * coeff);

        retMat.block(0,1) *= (1.0/M_alphaf);
        retMat.block(1,1) *= (1.0/M_alphaf);

        retMat += this->M_funProvider->getMass(time + dt*M_alphaf, solnpalphaf) * M_alpham;

        return retMat;
    });

    BV dersol = this->M_systemSolver.solve(fct, jac, initialGuess, status);

    if (status != 0)
        throw new Exception("Solver has not converged!");

    BV sol = computesolnp1(dersol, dt);

    std::vector<SolverStatistics> statistics = this->M_systemSolver.getSolverStatistics();
    this->dumpSolverStatistics(statistics, time+dt);

    return sol;
}

template <>
BlockVector<VectorEp>
GeneralizedAlphaMethod1stOrderPressure<VectorEp, MatrixEp>::
computeDerivative(const BlockVector<VectorEp>& solnp1, double& dt)
{
    BlockVector<VectorEp> dersolnpalpham;
    dersolnpalpham.hardCopy(M_prevDerivative);
    dersolnpalpham += (solnp1 - M_prevSolution - (M_prevDerivative * dt)) * (1.0 / (M_gamma * dt));

    dersolnpalpham.block(1) = (solnp1.block(1) - M_prevSolution.block(1));
    dersolnpalpham.block(1) *= (1.0 / dt);

    return dersolnpalpham;
}

}
