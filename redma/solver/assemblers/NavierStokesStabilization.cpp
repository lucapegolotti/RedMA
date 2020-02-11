#include "NavierStokesStabilization.hpp"

namespace RedMA
{

NavierStokesStabilization::
NavierStokesStabilization(const DataContainer& data,
                          SHP(FESPACE) fespaceVelocity,
                          SHP(FESPACE) fespacePressure,
                          SHP(ETFESPACE3) etfespaceVelocity,
                          SHP(ETFESPACE1) etfespacePressure) :
  M_timeOrder(data("time_discretization/order", 2)),
  M_dt(data("time_discretization/dt", 0.01)),
  M_velocityFESpace(fespaceVelocity),
  M_pressureFESpace(fespacePressure),
  M_velocityFESpaceETA(etfespaceVelocity),
  M_pressureFESpaceETA(etfespacePressure)
{
    std::string velocityOrder = data("fluid/velocity_order", "P1");
    if (!std::strcmp(velocityOrder.c_str(),"P1"))
        M_C_I = 30;
    else if (!std::strcmp(velocityOrder.c_str(),"P2"))
        M_C_I = 60;
    else if (!std::strcmp(velocityOrder.c_str(),"P3"))
        M_C_I = 120;
    else if (!std::strcmp(velocityOrder.c_str(),"P4"))
        M_C_I = 240;
    else
    {
        std::string msg = "Please implement a suitable value for ";
        msg += " M_C_I for your velocity FE order";
        throw Exception(msg);
    }
}

void
NavierStokesStabilization::
setDensityAndViscosity(const double& density, const double& viscosity)
{
    M_density = density;
    M_viscosity = viscosity;
}

}  // namespace RedMA
