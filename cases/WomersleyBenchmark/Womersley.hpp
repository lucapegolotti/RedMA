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

#ifndef WOMERSLEY_HPP
#define WOMERSLEY_HPP 1

#include <complex>

#include <lifev/core/LifeV.hpp>
#include <lifev/core/filter/GetPot.hpp>

namespace LifeV
{

class Womersley
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
    static double uexact_dtdt(const double& t, const double& x, const double& y,
                            const double& z, const ID& i);
    static double pexact (const double& t, const double& x, const double& y,
                        const double& z, const ID& i);

    static double pexact_dt(const double& t, const double& x, const double& y,
                          const double& z, const ID& i);

    // Initial velocity
    static double x0(const double& t, const double& x, const double& y,
                   const double& z, const ID& i);

    static double u0(const double& t, const double& x, const double& y,
                   const double& z, const ID& i);

    static double p0(const double& t, const double& x, const double& y,
                   const double& z, const ID& i);

    static double grad_u(const UInt& icoor, const double& t, const double& x, const double& y,
                       const double& z, const ID& i);

    static double grad_u_dt(const UInt& icoor, const double& t, const double& x, const double& y,
                            const double& z, const ID& i);

    static double fNeumann(const double& t, const double& x, const double& y,
                         const double& z, const ID& i);

    static double fNeumann_dt(const double& t, const double& x, const double& y,
                              const double& z, const ID& i);

    static double normalVector(const double& t, const double& x, const double& y,
                             const double& z, const ID& i);

    static double fShearStress(const double& t, const double& x, const double& y,
                             const double& z, const ID& i);

    static double fWallShearStress(const double& t, const double& x, const double& y,
                                 const double& z, const ID& i);
    static void setParamsFromGetPot(const GetPot& dataFile);
    static void computeQuantities();
    static void setPeriod(const double& newPeriod);
    static double getPeriod();
    static void showMe();

private:

    static Int  S_flagStrain;
    static double S_mu;
    static double S_nu;
    static double S_D;
    static double S_T;
    static double S_rho;
    static double S_W0;
    static double S_L;
    static double S_A;
    static double S_w;
    static std::complex<double> S_cj1, S_cy0, S_cy1, S_cj0p, S_cj1p, S_cy0p;
    static std::complex<double> S_cy1p , S_z1, S_b1 , S_ii, S_wi;
}; // class Womersley

} // namespace LifeV

#endif /* WOMERSLEY_HPP */
