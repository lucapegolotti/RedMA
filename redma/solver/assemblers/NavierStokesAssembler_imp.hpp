namespace RedMA
{

template <class InVectorType, class InMatrixType>
NavierStokesAssembler<InVectorType, InMatrixType>::
NavierStokesAssembler(const DataContainer& data, SHP(TreeNode) treeNode) :
  StokesAssembler<InVectorType, InMatrixType>(data, treeNode)
{
    this->M_name = "NavierStokesAssembler";
    M_useStabilization = this->M_data("assembler/use_stabilization", false);
}

template <class InVectorType, class InMatrixType>
void
NavierStokesAssembler<InVectorType, InMatrixType>::
setup()
{
    StokesAssembler<InVectorType, InMatrixType>::setup();

    if (M_useStabilization)
    {
        M_stabilization.reset(new SUPGStabilization(this->M_data,
                                                    this->M_velocityFESpace,
                                                    this->M_pressureFESpace,
                                                    this->M_velocityFESpaceETA,
                                                    this->M_pressureFESpaceETA));
        M_stabilization->setDensityAndViscosity(this->M_density, this->M_viscosity);
    }

}

template <class InVectorType, class InMatrixType>
BlockMatrix<InMatrixType>
NavierStokesAssembler<InVectorType, InMatrixType>::
getMass(const double& time, const BlockVector<InVectorType>& sol)
{
    BlockMatrix<InMatrixType> retMat;
    retMat.hardCopy(this->M_mass);
    if (M_useStabilization)
    {
        retMat += M_stabilization->getMass(sol, this->getForcingTerm(time));
        this->M_bcManager->apply0DirichletMatrix(retMat, this->getFESpaceBCs(),
                                                 this->getComponentBCs(), 1.0);
    }

    return retMat;
}

template <class InVectorType, class InMatrixType>
BlockMatrix<InMatrixType>
NavierStokesAssembler<InVectorType, InMatrixType>::
getMassJacobian(const double& time, const BlockVector<InVectorType>& sol)
{
    BlockMatrix<InMatrixType> retMat(this->M_nComponents, this->M_nComponents);
    if (M_useStabilization)
    {
        retMat += M_stabilization->getMassJac(sol, this->getForcingTerm(time));
        // we do it here because matrices from stabilization have no bcs
        this->M_bcManager->apply0DirichletMatrix(retMat, this->getFESpaceBCs(),
                                                 this->getComponentBCs(), 0.0);
    }

    return retMat;
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

    this->addNeumannBCs(retVec, time, sol);

    if (M_useStabilization)
        retVec -= M_stabilization->getResidual(sol, this->getForcingTerm(time));

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

    this->addConvectiveTermJacobianRightHandSide(sol, this->getZeroVector(), retMat);

    if (M_useStabilization)
    {
        retMat -= M_stabilization->getJac(sol, this->getForcingTerm(time));
        this->M_bcManager->apply0DirichletMatrix(retMat, this->getFESpaceBCs(),
                                                 this->getComponentBCs(), 0.0);
    }
    return retMat;
}


}
