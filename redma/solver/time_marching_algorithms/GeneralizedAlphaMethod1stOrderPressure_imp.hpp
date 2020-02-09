namespace RedMA
{

template <class InVectorType, class InMatrixType>
GeneralizedAlphaMethod1stOrderPressure<InVectorType, InMatrixType>::
GeneralizedAlphaMethod1stOrderPressure(const DataContainer& data) :
  aTimeMarchingAlgorithm<InVectorType, InMatrixType>(data),
  M_order(2)
{
}

template <class InVectorType, class InMatrixType>
GeneralizedAlphaMethod1stOrderPressure<InVectorType, InMatrixType>::
GeneralizedAlphaMethod1stOrderPressure(const DataContainer& data, SHP(FunProvider) funProvider) :
  aTimeMarchingAlgorithm<InVectorType, InMatrixType>(data, funProvider),
  M_order(2)
{
  setup(this->M_funProvider->getZeroVector());
}

template <class InVectorType, class InMatrixType>
void
GeneralizedAlphaMethod1stOrderPressure<InVectorType, InMatrixType>::
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
GeneralizedAlphaMethod1stOrderPressure<InVectorType, InMatrixType>::
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
GeneralizedAlphaMethod1stOrderPressure<InVectorType, InMatrixType>::
computedersolnpalpham(BlockVector<InVectorType> dersol)
{
    BlockVector<InVectorType> dersolnpalpham;
    dersolnpalpham.hardCopy(M_prevDerivative);
    dersolnpalpham += (dersol - M_prevDerivative) * M_alpham;

    return dersolnpalpham;
}

template <class InVectorType, class InMatrixType>
void
GeneralizedAlphaMethod1stOrderPressure<InVectorType, InMatrixType>::
shiftSolutions(const BlockVector<InVectorType>& sol)
{
    double dt = this->M_data("time_discretization/dt", 0.01);
    M_prevDerivative = computeDerivative(sol, dt);
    M_prevSolution = sol;
}

}
