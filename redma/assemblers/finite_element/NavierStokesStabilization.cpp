#include "NavierStokesStabilization.hpp"

namespace RedMA
{

NavierStokesStabilization::
NavierStokesStabilization(const DataContainer& data,
                          shp<FESPACE> fespaceVelocity,
                          shp<FESPACE> fespacePressure,
                          shp<ETFESPACE3> etfespaceVelocity,
                          shp<ETFESPACE1> etfespacePressure,
                          EPETRACOMM comm) :
  M_velocityFESpace(fespaceVelocity),
  M_pressureFESpace(fespacePressure),
  M_velocityFESpaceETA(etfespaceVelocity),
  M_pressureFESpaceETA(etfespacePressure),
  M_comm(comm)
{
    M_velocityOrder = data("fluid/velocity_order", "P1");
    M_verbose = data.getVerbose();
}

void
NavierStokesStabilization::
setDensityAndViscosity(const double& density,
                       const double& viscosity)
{
    M_density = density;
    M_viscosity = viscosity;
}

}  // namespace RedMA
