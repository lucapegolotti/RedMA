namespace RedMA
{

template <class InVectorType, class InMatrixType>
NavierStokesAssembler<InVectorType, InMatrixType>::
NavierStokesAssembler(const DataContainer& data, SHP(TreeNode) treeNode) :
  StokesAssembler<InVectorType, InMatrixType>(data, treeNode),
  M_useStabilization(false),
  M_stabilization(nullptr)
{
    this->M_name = "NavierStokesAssembler";
}

template <class InVectorType, class InMatrixType>
void
NavierStokesAssembler<InVectorType, InMatrixType>::
setup()
{
    StokesAssembler<InVectorType, InMatrixType>::setup();

    std::string stabilizationType = this->M_data("assembler/stabilization", "none");

    if (!std::strcmp(stabilizationType.c_str(),"supg"))
    {
        M_stabilization.reset(new SUPGStabilization(this->M_data,
                                                    this->M_velocityFESpace,
                                                    this->M_pressureFESpace,
                                                    this->M_velocityFESpaceETA,
                                                    this->M_pressureFESpaceETA));
        M_stabilization->setDensityAndViscosity(this->M_density, this->M_viscosity);
    }
    else if (!std::strcmp(stabilizationType.c_str(),"vms"))
    {
        M_stabilization.reset(new VMSStabilization(this->M_data,
                                                   this->M_velocityFESpace,
                                                   this->M_pressureFESpace,
                                                   this->M_velocityFESpaceETA,
                                                   this->M_pressureFESpaceETA));
        M_stabilization->setDensityAndViscosity(this->M_density, this->M_viscosity);
    }
    else if (!std::strcmp(stabilizationType.c_str(),"hf"))
    {
        M_stabilization.reset(new HFStabilization(this->M_data,
                                                  this->M_velocityFESpace,
                                                  this->M_pressureFESpace,
                                                  this->M_velocityFESpaceETA,
                                                  this->M_pressureFESpaceETA));
        M_stabilization->setDensityAndViscosity(this->M_density, this->M_viscosity);
    }
    else if (!std::strcmp(stabilizationType.c_str(),"none"))
    {

    }
    else
        throw new Exception("Type of stabilization not implemented!");

    M_useStabilization = (M_stabilization != nullptr);
}

}
