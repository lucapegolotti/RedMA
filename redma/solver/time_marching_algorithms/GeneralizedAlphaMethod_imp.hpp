namespace RedMA
{

template <class InVectorType, class InMatrixType>
GeneralizedAlphaMethod<InVectorType, InMatrixType>::
GeneralizedAlphaMethod(const DataContainer& data) :
  aTimeMarchingAlgorithm<InVectorType, InMatrixType>(data),
  M_order(2)
{

}

template <class InVectorType, class InMatrixType>
GeneralizedAlphaMethod<InVectorType, InMatrixType>::
GeneralizedAlphaMethod(const DataContainer& data, SHP(FunProvider) funProvider) :
  aTimeMarchingAlgorithm<InVectorType, InMatrixType>(data, funProvider),
  M_order(2)
{
  setup(this->M_funProvider->getZeroVector());
}

template <class InVectorType, class InMatrixType>
void
GeneralizedAlphaMethod<InVectorType, InMatrixType>::
setup(const BlockVector<InVectorType>& zeroVector)
{
    M_prevSolution.hardCopy(zeroVector);
    M_prevDerivative.hardCopy(zeroVector);

    M_rhoinf = this->M_data("time_discretization/rhoinf", 0.5);

    M_alpham = 0.5 * (3.0 - M_rhoinf) / (1.0 + M_rhoinf);
    M_alphaf = 1.0 / (1.0 + M_rhoinf);
    M_gamma = 0.5 + M_alpham - M_alphaf;
}

template <class InVectorType, class InMatrixType>
BlockVector<InVectorType>
GeneralizedAlphaMethod<InVectorType, InMatrixType>::
advance(const double& time, double& dt, int& status)
{
    typedef BlockVector<InVectorType>               BV;
    typedef BlockMatrix<InMatrixType>               BM;

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
        retMat += this->M_funProvider->getMass(time + dt*M_alphaf, solnpalphaf) * M_alpham;
        retMat += (this->M_funProvider->getMassJacobian(time + dt*M_alphaf, solnpalphaf) * coeff);

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

template <class InVectorType, class InMatrixType>
BlockVector<InVectorType>
GeneralizedAlphaMethod<InVectorType, InMatrixType>::
computesolnp1(BlockVector<InVectorType> dersol, const double& dt)
{
    BlockVector<InVectorType> solnp1;
    solnp1.hardCopy(M_prevSolution);
    solnp1 += (M_prevDerivative * dt);
    solnp1 += (dersol - M_prevDerivative) * (M_gamma * dt);

    return solnp1;
}

template <class InVectorType, class InMatrixType>
BlockVector<InVectorType>
GeneralizedAlphaMethod<InVectorType, InMatrixType>::
computesolnpalphaf(BlockVector<InVectorType> solnp1)
{
    BlockVector<InVectorType> solnpalphaf;
    solnpalphaf.hardCopy(M_prevSolution);
    solnpalphaf += (solnp1 - M_prevSolution) * M_alphaf;

    return solnpalphaf;
}

template <class InVectorType, class InMatrixType>
BlockVector<InVectorType>
GeneralizedAlphaMethod<InVectorType, InMatrixType>::
computedersolnpalpham(BlockVector<InVectorType> dersol)
{
    BlockVector<InVectorType> dersolnpalpham;
    dersolnpalpham.hardCopy(M_prevDerivative);
    dersolnpalpham += (dersol - M_prevDerivative) * M_alpham;

    return dersolnpalpham;
}

template <class InVectorType, class InMatrixType>
BlockVector<InVectorType>
GeneralizedAlphaMethod<InVectorType, InMatrixType>::
computeDerivative(const BlockVector<InVectorType>& solnp1, double& dt)
{
    BlockVector<InVectorType> dersolnpalpham;
    dersolnpalpham.hardCopy(M_prevDerivative);
    dersolnpalpham += (solnp1 - M_prevSolution - (M_prevDerivative * dt)) * (1.0 / (M_gamma * dt));
    return dersolnpalpham;
}

template <class InVectorType, class InMatrixType>
void
GeneralizedAlphaMethod<InVectorType, InMatrixType>::
shiftSolutions(const BlockVector<InVectorType>& sol)
{
    double dt = this->M_data("time_discretization/dt", 0.01);
    M_prevDerivative = computeDerivative(sol, dt);
    M_prevSolution = sol;
}

}
