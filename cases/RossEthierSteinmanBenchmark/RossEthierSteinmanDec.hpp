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

#ifndef ROSS_ETHIER_STEINMAN_DEC_HPP
#define ROSS_ETHIER_STEINMAN_DEC_HPP 1

#include <lifev/core/LifeV.hpp>
#include <lifev/core/filter/GetPot.hpp>

namespace LifeV
{

class RossEthierSteinmanDec
{
public:
    static double f(const double& t, const double& x, const double& y,
                    const double& z, const ID& i);

    static double xexact(const double& t, const double& x, const double& y,
                         const double& z, const ID& i);

    static double uexact(const double& t, const double& x, const double& y,
                         const double& z, const ID& i);

    static double uexact_dt(const double& t, const double& x, const double& y,
                            const double& z, const ID& i);

    static double uderexact(const double& t, const double& x, const double& y,
                            const double& z, const ID& i);

    static double pexact(const double& t, const double& x, const double& y,
                         const double& z, const ID& i);

    static double pexact_dt(const double& t, const double& x, const double& y,
                            const double& z, const ID& i);

    static double grad_u(const UInt& icoor, const double& t, const double& x, const double& y,
                         const double& z, const ID& i);

    static double grad_u_dt(const UInt& icoor, const double& t, const double& x, const double& y,
                            const double& z, const ID& i);

    // Initial velocity
    static double x0(const double& t, const double& x, const double& y,
                     const double& z, const ID& i);

    static double fNeumann(const double& t, const double& x, const double& y,
                           const double& z, const ID& i);

    static double fNeumann_dt(const double& t, const double& x, const double& y,
                              const double& z, const ID& i);

    static double normalVector(const double& t, const double& x, const double& y,
                               const double& z, const ID& i);

    static void setParamsFromGetPot(const GetPot& dataFile);
    static void setA (const double& aValue);
    static void setD (const double& dValue);
    static void setViscosity (const double& mu);
    static void setDensity (const double& rho);
    static void setFlagStrain (const Int& flagValue);

private:

    static double S_a;
    static double S_d;
    static double S_L;
    static double S_mu;
    static double S_rho;
    static double S_nu;
    static Int  S_flagStrain;
}; // class RossEthierSteinmanDec

} // namespace LifeV

#endif /* ROSS_ETHIER_STEINMAN_DEC_HPP */
