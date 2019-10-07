// Reduced Modeling of Arteries (RedMA)
// Copyright (C) 2019  Luca Pegolotti
//
// RedMA is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// RedMA is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with this program.  If not, see <http://www.gnu.org/licenses/>.

#ifndef UDFNS_HPP
#define UDFNS_HPP

// LifeV includes
#include <lifev/core/LifeV.hpp>

#include "Womersley.hpp"

namespace LifeV
{

double fZero(const double& /*t*/, const double& /*x*/, const double& /*y*/,
             const double& /*z*/, const ID& /*i*/)
{
	return 0.0;
}

double fOne(const double& /*t*/, const double& /*x*/, const double& /*y*/,
            const double& /*z*/, const ID& /*i*/)
{
	return 1.0;
}

double pressureWomerTime(const double& t, const double& /*x*/, const double& y,
                         const double& z, const ID& /*i*/)
{
  return  -1.*Womersley::pexact(t, 0.0, y, z, 0);
}

double pressuredtWomerTime(const double& t, const double& /*x*/, const double& y,
                           const double& z, const ID& /*i*/)
{
    return  -1.*Womersley::pexact_dt(t, 0.0, y, z, 0);
}

double velocityWomerTime(const double& t, const double& /*x*/, const double& y,
                         const double& z, const ID& i)
{
   return Womersley::uexact(t, 0.0, y, z, i);
}

double velocitydtWomerTime(const double& t, const double& /*x*/, const double& y,
                           const double& z, const ID& i)
{
    return Womersley::uexact_dt(t, 0.0, y, z, i);
}

double velocitydtdtWomerTime(const double& t, const double& /*x*/, const double& y,
                             const double& z, const ID& i)
{
    return Womersley::uexact_dtdt(t, 0.0, y, z, i);
}

} // namespace LifeV

#endif
