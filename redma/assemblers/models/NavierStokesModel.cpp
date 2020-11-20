#include "NavierStokesModel.hpp"

namespace RedMA
{

NavierStokesModel::
NavierStokesModel(const DataContainer& data, shp<TreeNode> treeNode)
  // StokesModel(data, treeNode),
  // M_useStabilization(false)
  // M_stabilization(nullptr)
{
    // this->M_name = "NavierStokesModel";
}

// void
// NavierStokesAssembler::
// setup()
// {
//     StokesAssembler::setup();
//
//     // std::string stabilizationType = this->M_data("assembler/stabilization", "none");
//     //
//     // if (!std::strcmp(stabilizationType.c_str(),"supg"))
//     // {
//     //     M_stabilization.reset(new SUPGStabilization(this->M_data,
//     //                                                 this->M_velocityFESpace,
//     //                                                 this->M_pressureFESpace,
//     //                                                 this->M_velocityFESpaceETA,
//     //                                                 this->M_pressureFESpaceETA));
//     //     M_stabilization->setDensityAndViscosity(this->M_density, this->M_viscosity);
//     // }
//     // else if (!std::strcmp(stabilizationType.c_str(),"vms"))
//     // {
//     //     M_stabilization.reset(new VMSStabilization(this->M_data,
//     //                                                this->M_velocityFESpace,
//     //                                                this->M_pressureFESpace,
//     //                                                this->M_velocityFESpaceETA,
//     //                                                this->M_pressureFESpaceETA));
//     //     M_stabilization->setDensityAndViscosity(this->M_density, this->M_viscosity);
//     // }
//     // else if (!std::strcmp(stabilizationType.c_str(),"hf"))
//     // {
//     //     M_stabilization.reset(new HFStabilization(this->M_data,
//     //                                               this->M_velocityFESpace,
//     //                                               this->M_pressureFESpace,
//     //                                               this->M_velocityFESpaceETA,
//     //                                               this->M_pressureFESpaceETA));
//     //     M_stabilization->setDensityAndViscosity(this->M_density, this->M_viscosity);
//     // }
//     // else if (!std::strcmp(stabilizationType.c_str(),"none"))
//     // {
//     //
//     // }
//     // else
//     //     throw new Exception("Type of stabilization not implemented!");
//
//     M_useStabilization = (M_stabilization != nullptr);
// }

}
