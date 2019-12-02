#include "NavierStokesParameterHandler.hpp"

namespace RedMA
{

NavierStokesParameterHandler::
NavierStokesParameterHandler(const param_Type& _muMin,
                             const param_Type& _muMax,
                             bool _verbose) :
  ParameterHandler(_muMin, _muMax, _verbose)
{
	setTheta();
}

void
NavierStokesParameterHandler::
setTheta()
{
}


} // namespace RedMA
