#include <lifev/core/LifeV.hpp>
#include "Womersley.hpp"
#include <lifev/navier_stokes_blocks/function/bessel/bessel.hpp>

namespace LifeV
{

double
Womersley::
uexact(const double& t, const double& /*x*/, const double& y,
       const double& z, const ID& i)
{
    double r = std::sqrt (z * z + y * y);
    std::complex<double> z2, b2;
    z2 = 2.*r / S_D * S_z1;
    bessel::cbessjy01 (z2, b2, S_cj1, S_cy0, S_cy1, S_cj0p, S_cj1p, S_cy0p, S_cy1p);
    double u = real(S_A / S_L / S_rho / S_wi * (1. - b2 / S_b1) * std::exp (S_wi * t) );
    switch (i)
    {
        case 0:  //u_1
            return u;//-4*x*y*y; //u-4*x*y*y;
        case 1:  //u_2
            return 0.;//y*y*y;
        case 2:
            return 0.;
        default:
            exit (1);
    }
    return 1.;
}

double
Womersley::
uexact_dt(const double& t, const double& /*x*/, const double& y,
          const double& z, const ID& i)
{
    double r = std::sqrt (z * z + y * y);
    std::complex<double> z2, b2;
    z2 = 2.*r / S_D * S_z1;
    bessel::cbessjy01 (z2, b2, S_cj1, S_cy0, S_cy1, S_cj0p, S_cj1p, S_cy0p, S_cy1p);
    double ut = real (S_A / S_L / S_rho * (1. - b2 / S_b1) * std::exp (S_wi * t) );
    switch (i)
    {
        case 0:  //u_1
            return ut;//-4*x*y*y; //u-4*x*y*y;
        case 1:  //u_2
            return 0.;//y*y*y;
        case 2:
            return 0.;
        default:
            exit (1);
    }
    return 1.;
}

double
Womersley::
uexact_dtdt(const double& t, const double& /*x*/,
            const double& y, const double& z, const ID& i)
{
    double r = std::sqrt (z * z + y * y);
    std::complex<double> z2, b2;
    z2 = 2.*r / S_D * S_z1;
    bessel::cbessjy01 (z2, b2, S_cj1, S_cy0, S_cy1, S_cj0p, S_cj1p, S_cy0p, S_cy1p);
    double utt =  real(S_wi * S_A / S_L / S_rho * (1. - b2 / S_b1) * std::exp (S_wi * t));
    switch (i)
    {
        case 0:  //u_1
            return utt;//-4*x*y*y; //u-4*x*y*y;
        case 1:  //u_2
            return 0.;//y*y*y;
        case 2:
            return 0.;
        default:
            exit (1);
    }
    return 1.;
}

double
Womersley::
pexact(const double& t, const double& x,
       const double& /*y*/, const double& /*z*/, const ID& /* i */)
{
    return S_A / S_L * (S_L - x) * std::cos (S_w * t);
}

double
Womersley::
pexact_dt(const double& t, const double& x,
          const double& /*y*/, const double& /*z*/, const ID& /* i */)
{
    return S_A / S_L * (S_L - x) * S_w * std::cos (S_w * t);
}

double
Womersley::
grad_u(const UInt& icoor, const double& t, const double& /*x*/,
       const double& y, const double& z, const ID& i )
{
    double r = std::sqrt (y * y + z * z);
    std::complex<double> z2, b2;
    z2 = 2.*r / S_D * S_z1;
    bessel::cbessjy01(z2, b2, S_cj1, S_cy0, S_cy1, S_cj0p, S_cj1p, S_cy0p, S_cy1p);
    b2 = -2. / S_D * S_z1 * S_cj0p;
    double u_r = real(S_A / S_L / S_rho / S_wi * +b2 / S_b1 * std::exp (S_wi * t) );

    switch (icoor)
    {
        case 0:
            switch (i)
            {
                case 0:
                    return 0.;
                case 1:
                    return 0.;
                case 2:
                    return 0.;
                default:
                    exit (1);
            }
        case 1:   // u_y
            switch (i)
            {
                case 0:
                    return u_r / r * y;
                case 1:
                    return 0.;
                case 2:
                    return 0.;
                default:
                    exit (1);
            }
        case 2:
            switch (i)
            {
                case 0:
                    return u_r / r * z;
                case 1:
                    return 0.;
                case 2:
                    return 0.;
                default:
                    exit (1);
            }
        default:
            exit (1);
    }
    return 1.;
}

double
Womersley::
f(const double& /* t */, const double&  /*x*/, const double&  /*y*/,
  const double& /* z */, const ID& i )
{
    switch (i)
    {
        case 0:
            return 0.;
        case 1:
            return 0.;
        case 2:
            return 0.;
        default:
            exit (1);
    }
    return 1.;
}

double
Womersley::
xexact(const double& t, const double& x, const double& y,
       const double& z, const ID& i)
{
    switch (i)
    {
        case 0:  //u_1
        case 1:  //u_2
            return uexact (t, x, y, z, i);
        case 2:  //pressure
            return pexact (t, x, y, z, 1);
            break;
        default:
            exit (1);
    }
    return 1.;
}

double
Womersley::
x0(const double& t, const double& x,
   const double& y, const double& z, const ID& i)
{
    return xexact (t, x, y, z, i);
}

// Initial velocity
double
Womersley::
u0(const double& t, const double& x, const double& y,
   const double& z, const ID& i)
{
    return uexact (t, x, y, z, i);
}

// Initial pressure
double
Womersley::
p0(const double& t, const double& x, const double& y,
   const double& z, const ID& /*i*/ )
{
    return pexact (t, x, y, z, 0);
}

//we suppose that the problem geometry is the cylinder having axis x, origin (0,0,0), diameter D and height L
double
Womersley::
fNeumann(const double& t, const double& x, const double& y,
         const double& z, const ID& i )
{
    double r = std::sqrt (y * y + z * z);
    double n[3] = {0., 0., 0.};
    double out = 0.;
    if        (x < 1e-6 / S_L )
    {
        n[0] = -1.;
    }
    else if (x >  S_L * (1 - 1e-6) )
    {
        n[0] =  1.;
    }
    else if (r > S_D / 2.*0.9 )
    {
        n[1] = y / r;
        n[2] = z / r;
    }
    else
    {
        // std::cout << "strange point: x=" << x << " y=" << y << " z=" << z
        //          << std::endl;
    }

    for (UInt k = 0; k < nDimensions; k++) //mu gradu n
    {
        out += S_mu * grad_u (k, t, x, y, z, i) * n[k];
    }

    if (S_flagStrain)
        for (UInt k = 0; k < nDimensions; k++) //mu gradu^T n
        {
            out += S_mu * grad_u (i, t, x, y, z, k) * n[k];
        }

    out -= pexact (t, x, y, z, i) * n[i];

    return  out;
}

double
Womersley::
normalVector(const double& /*t*/, const double& x, const double& y,
             const double& z, const ID& i )
{
    double r = std::sqrt (y * y + z * z);
    double n[3] = {0., 0., 0.};
    if        (x < 1e-6 / S_L )
    {
        n[0] = -1.;
    }
    else if (x >  S_L * (1 - 1e-6) )
    {
        n[0] =  1.;
    }
    else if (r > S_D / 2.*0.9 )
    {
        n[1] = y / r;
        n[2] = z / r;
    }
    else
    {
        // std::cout << "strange point: x=" << x << " y=" << y << " z=" << z
        //          << std::endl;
    }
    return n[i];
}

//we suppose that the problem geometry is the cylinder having axis x, origin (0,0,0), diameter D and height L
double
Womersley::
fShearStress(const double& t, const double& x, const double& y,
             const double& z, const ID& i)
{
    double r = std::sqrt (y * y + z * z);
    double n[3] = {0., 0., 0.};
    double out = 0.;
    if        (x < 1e-6 / S_L )
    {
        n[0] = -1.;
    }
    else if (x >  S_L * (1 - 1e-6) )
    {
        n[0] =  1.;
    }
    else if (r > S_D / 2.*0.9 )
    {
        n[1] = y / r;
        n[2] = z / r;
    }
    else
    {
        // std::cout << "strange point: x=" << x << " y=" << y << " z=" << z
        //          << std::endl;
    }

    for (UInt k = 0; k < nDimensions; k++) //mu gradu n
    {
        out += S_mu * grad_u (k, t, x, y, z, i) * n[k];
    }

    if (S_flagStrain)
        for (UInt k = 0; k < nDimensions; k++) //mu gradu^T n
        {
            out += S_mu * grad_u (i, t, x, y, z, k) * n[k];
        }

    return  out;
}

//we suppose that the problem geometry is the cylinder having axis x, origin (0,0,0), diameter D and height L
double
Womersley::
fWallShearStress(const double& t, const double& x, const double& y,
                 const double& z, const ID& i )
{
    double r = std::sqrt (y * y + z * z);
    double n[3] = {0., 0., 0.};
    if        (x < 1e-6 / S_L )
    {
        n[0] = -1.;
    }
    else if (x >  S_L * (1 - 1e-6) )
    {
        n[0] =  1.;
    }
    else if (r > S_D / 2.*0.9 )
    {
        n[1] = y / r;
        n[2] = z / r;
    }
    else
    {
        //std::cout << "strange point: x=" << x << " y=" << y << " z=" << z
        //          << std::endl;
    }

    // wss = s_n - (s_n \cdot n ) n
    double wss = 0.;
    double s_n_n (0.);
    // this is the i-component of the normal stress
    wss = fShearStress (t, x, y, z, i);

    for (UInt k = 0; k < nDimensions; k++) // (s_n \cdot n )
    {
        s_n_n += fShearStress (t, x, y, z, k) * n[k];
    }

    wss -= s_n_n * n[i];

    return  wss;
}

void
Womersley::
setParamsFromGetPot(const GetPot& dataFile )
{
    S_mu = dataFile("fluid/viscosity", 1.);
    S_flagStrain = dataFile("fluid/flag_strain", 0);
    S_nu = S_mu / dataFile("fluid/density", 1.);
    S_D = dataFile("fluid/D", 6.);
    S_T = dataFile("fluid/T", 15.);
    S_rho = dataFile("fluid/density", 1.);
    S_L = dataFile("fluid/L", 1.);
    S_A = dataFile("fluid/A", 16 * S_mu / S_D / S_D * S_L * S_L);
    computeQuantities();
}

void
Womersley::
computeQuantities()
{
    S_W0 = S_D / 2 * std::sqrt (2 * M_PI / S_nu / S_T);
    S_ii = std::complex<double> (0., 1.);
    S_w = 2.*M_PI / S_T;
    S_wi = S_w * S_ii;
    S_z1 = S_W0 * std::pow (S_ii, 1.5);
    bessel::cbessjy01 (S_z1, S_b1, S_cj1, S_cy0, S_cy1, S_cj0p, S_cj1p, S_cy0p, S_cy1p);
}

double
Womersley::
getPeriod()
{
    return S_T;
}

void
Womersley::
setPeriod(const double& newPeriod)
{
    S_T = newPeriod;
    computeQuantities();
}

void
Womersley::
showMe()
{
    std::cout << "Kynetic viscosity " << S_nu << std::endl;
    std::cout << "M_PIpe radius " << S_D / 2. << " M_PIpe lenght " << S_L << std::endl;
    std::cout << "Oscillation period " << S_T << std::endl;
    std::cout << "Pressure Drop " << S_A << std::endl;
    std::cout << "Womersley Number " << S_D / 2.*std::sqrt (S_w / S_nu) << std::endl;
}

double Womersley::S_nu;
double Womersley::S_mu;
Int  Womersley::S_flagStrain;
double Womersley::S_D;
double Womersley::S_rho;
double Womersley::S_T;
double Womersley::S_W0;
double Womersley::S_L;
double Womersley::S_A;
double Womersley::S_w;
std::complex<double> Womersley::S_cj1, Womersley::S_cy0, Womersley::S_cy1,
        Womersley::S_cj0p, Womersley::S_cj1p, Womersley::S_cy0p, Womersley::S_cy1p ,
        Womersley::S_z1, Womersley::S_b1 , Womersley::S_ii, Womersley::S_wi;
} // namespace LifeV
