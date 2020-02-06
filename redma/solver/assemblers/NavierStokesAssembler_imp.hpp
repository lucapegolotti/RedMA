namespace RedMA
{

template <class InVectorType, class InMatrixType>
NavierStokesAssembler<InVectorType, InMatrixType>::
NavierStokesAssembler(const DataContainer& data, SHP(TreeNode) treeNode) :
  StokesAssembler<InVectorType, InMatrixType>(data, treeNode)
{

}

template <class InVectorType, class InMatrixType>
BlockVector<InVectorType>
NavierStokesAssembler<InVectorType, InMatrixType>::
getRightHandSide(const double& time, const BlockVector<InVectorType>& sol)
{
    BlockVector<InVectorType> retVec;
    BlockMatrix<InMatrixType> systemMatrix;

    systemMatrix.resize(this->M_nComponents, this->M_nComponents);
    systemMatrix += this->M_stiffness;
    systemMatrix += this->M_divergence;
    systemMatrix *= (-1.0);

    this->addConvectiveMatrixRightHandSide(sol, systemMatrix);

    retVec.softCopy(systemMatrix * sol);

    // treatment of Dirichlet bcs if needed
    bool useLifting = this->M_data("bc_conditions/lifting", true);

    if (useLifting)
    {
        BlockVector<InVectorType> lifting = this->computeLifting(time);
        retVec += (systemMatrix * lifting);

        double dt = this->M_data("time_discretization/dt", 0.01);
        BlockVector<InVectorType> liftingDt = this->M_TMAlifting->computeDerivative(lifting, dt);
        retVec -= (this->M_mass * liftingDt);
    }

    this->addNeumannBCs(retVec, time, sol);

    this->apply0DirichletBCs(retVec);
    return retVec;
}

template <class InVectorType, class InMatrixType>
BlockMatrix<InMatrixType>
NavierStokesAssembler<InVectorType, InMatrixType>::
getJacobianRightHandSide(const double& time,
                         const BlockVector<InVectorType>& sol)
{
    BlockMatrix<InMatrixType> retMat;
    retMat = StokesAssembler<InVectorType,InMatrixType>::getJacobianRightHandSide(time, sol);

    BlockVector<InVectorType> lifting;

    bool useLifting = this->M_data("bc_conditions/lifting", true);

    if (useLifting)
        lifting = this->computeLifting(time);


    this->addConvectiveTermJacobianRightHandSide(sol, lifting, retMat);
    return retMat;
}


}
