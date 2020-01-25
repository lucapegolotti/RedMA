namespace RedMA
{

template <class InVectorType, class InMatrixType>
StokesAssembler<InVectorType, InMatrixType>::
StokesAssembler(const GetPot& datafile,
                SHP(BuildingBlock) buildingBlock) :
  aAssembler<InVectorType, InMatrixType>(datafile),
  M_comm(buildingBlock->getComm()),
  M_buildingBlock(buildingBlock),
  M_nComponents(2)
{
    M_density = this->M_datafile("fluid/density", 1.0);
    M_viscosity = this->M_datafile("fluid/viscosity", 0.035);

    // we check if building block is inlet. If so, we increase the number of
    // components of the system in order to impose the boundary conditions weakly
    if (!M_buildingBlock->getIsChild())
    {
        M_nComponents = 3;
    }
}

template <class InVectorType, class InMatrixType>
void
StokesAssembler<InVectorType,InMatrixType>::
setup()
{
    initializeFEspaces();

    assembleStiffness();
    assembleMass();
}

template <class InVectorType, class InMatrixType>
void
StokesAssembler<InVectorType, InMatrixType>::
initializeFEspaces()
{
    // initialize fespace velocity
    std::string orderVelocity = this->M_datafile("fluid/velocity_order", "P2");

    M_velocityFESpace.reset(new FESPACE(M_buildingBlock->getMesh(),
                                        orderVelocity, 3, M_comm));
    M_velocityFESpaceETA.reset(new ETFESPACE3(M_velocityFESpace->mesh(),
                                            &(M_velocityFESpace->refFE()),
                                              M_comm));

    // initialize fespace velocity
    std::string orderPressure = this->M_datafile("fluid/pressure_order", "P1");

    M_pressureFESpace.reset(new FESPACE(M_buildingBlock->getMesh(),
                                        orderPressure, 1, M_comm));
    M_pressureFESpaceETA.reset(new ETFESPACE1(M_pressureFESpace->mesh(),
                                            &(M_pressureFESpace->refFE()),
                                              M_comm));

}

template <class InVectorType, class InMatrixType>
void
StokesAssembler<InVectorType, InMatrixType>::
exportSolution(const double& t)
{

}

template <class InVectorType, class InMatrixType>
void
StokesAssembler<InVectorType, InMatrixType>::
postProcess()
{

}

template <class InVectorType, class InMatrixType>
BlockMatrix<InMatrixType>
StokesAssembler<InVectorType, InMatrixType>::
getMass(const double& time, const BlockVector<InVectorType>& sol)
{
    return M_mass;
}

template <class InVectorType, class InMatrixType>
BlockVector<InVectorType>
StokesAssembler<InVectorType, InMatrixType>::
getRightHandSide(const double& time, const BlockVector<InVectorType>& sol)
{
    BlockVector<InVectorType> retVec;
    // build right hand side here ..

    return retVec;
}

template <class InVectorType, class InMatrixType>
BlockMatrix<InMatrixType>
StokesAssembler<InVectorType, InMatrixType>::
getJacobianRightHandSide(const double& time, const BlockVector<InVectorType>& sol)
{
    BlockMatrix<InMatrixType> retMat;

    return retMat;
}

}
