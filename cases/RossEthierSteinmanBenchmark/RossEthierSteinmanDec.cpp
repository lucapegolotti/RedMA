#if 1

#include "RossEthierSteinmanDec.hpp"

namespace LifeV
{

double
RossEthierSteinmanDec::
uexact(const double& t, const double& x, const double& y,
       const double& z, const ID& i)
{
    double e = std::exp(-S_d * S_d * S_nu * t);

    switch(i)
    {
        case 0:
            return -S_a * e *(std::exp(S_a * x) * std::sin(S_a * y + S_d * z) + std::exp(S_a * z) * std::cos(S_a * x + S_d * y));
        case 1:
            return -S_a * e *(std::exp(S_a * y) * std::sin(S_a * z + S_d * x) + std::exp(S_a * x) * std::cos(S_a * y + S_d * z));
        case 2:
            return -S_a * e *(std::exp(S_a * z) * std::sin(S_a * x + S_d * y) + std::exp(S_a * y) * std::cos(S_a * z + S_d * x));
        default:
            exit(1);
    }
    return 1.;
}

double
RossEthierSteinmanDec::
uexact_dt(const double& t, const double& x, const double& y,
          const double& z, const ID& i)
{
    double e = -std::exp(-S_d * S_d * S_nu * t) * S_d * S_d * S_nu;

    switch(i)
    {
        case 0:
            return -S_a * e *(std::exp(S_a * x) * std::sin(S_a * y + S_d * z) + std::exp(S_a * z) * std::cos(S_a * x + S_d * y));
        case 1:
            return -S_a * e *(std::exp(S_a * y) * std::sin(S_a * z + S_d * x) + std::exp(S_a * x) * std::cos(S_a * y + S_d * z));
        case 2:
            return -S_a * e *(std::exp(S_a * z) * std::sin(S_a * x + S_d * y) + std::exp(S_a * y) * std::cos(S_a * z + S_d * x));
        default:
            exit(1);
    }
    return 1.;
}

double
RossEthierSteinmanDec::
uderexact(const double& t, const double& x, const double& y,
          const double& z, const ID& i)
{

    if(i < 3)
    {
        return -S_d * S_d * S_nu * xexact(t, x, y, z, i);
    }
    else
    {
        return 0.;
    }
}

double
RossEthierSteinmanDec::
pexact(const double& t, const double& x, const double& y,
                                             const double& z,
                                             const ID& /* i */)
{
    return - S_rho * S_a * S_a / 2. * std::exp(-2.*S_d * S_d * S_nu * t) *
         (std::exp(2.*S_a * x) + std::exp(2.*S_a * y) + std::exp(2.*S_a * z) +
             2. * std::sin(S_a * x + S_d * y) * std::cos(S_a * z + S_d * x) * std::exp(S_a *(y + z)) +
             2. * std::sin(S_a * y + S_d * z) * std::cos(S_a * x + S_d * y) * std::exp(S_a *(z + x)) +
             2. * std::sin(S_a * z + S_d * x) * std::cos(S_a * y + S_d * z) * std::exp(S_a *(x + y)));
}

double
RossEthierSteinmanDec::
pexact_dt(const double& t, const double& x, const double& y,
                                             const double& z,
                                             const ID& /* i */)
{
    return - S_rho * S_a * S_a / 2. * std::exp(-2.*S_d * S_d * S_nu * t) * (-2.) * S_d * S_d * S_nu *
         (std::exp(2.*S_a * x) + std::exp(2.*S_a * y) + std::exp(2.*S_a * z) +
             2. * std::sin(S_a * x + S_d * y) * std::cos(S_a * z + S_d * x) * std::exp(S_a *(y + z)) +
             2. * std::sin(S_a * y + S_d * z) * std::cos(S_a * x + S_d * y) * std::exp(S_a *(z + x)) +
             2. * std::sin(S_a * z + S_d * x) * std::cos(S_a * y + S_d * z) * std::exp(S_a *(x + y)));
}

double
RossEthierSteinmanDec::
grad_u(const UInt& icoor, const double& t, const double& x, const double& y, const double& z, const ID& i)
{
    double e = std::exp(-S_d * S_d * S_nu * t);
    switch(icoor)
    {
        case 0:    // u_x
            switch(i)
            {
                case 0:
                    return -S_a * e *(S_a * std::exp(S_a * x) * std::sin(S_a * y + S_d * z) - S_a * std::exp(S_a * z) * std::sin(S_a * x + S_d * y));
                case 1:
                    return -S_a * e *(S_d * std::exp(S_a * y) * std::cos(S_a * z + S_d * x) + S_a * std::exp(S_a * x) * std::cos(S_a * y + S_d * z));
                case 2:
                    return -S_a * e *(S_a * std::exp(S_a * z) * std::cos(S_a * x + S_d * y) - S_d * std::exp(S_a * y) * std::sin(S_a * z + S_d * x));
                default:
                    exit(1);
            }
        case 1:   // u_y
            switch(i)
            {
                case 0:
                    return -S_a * e *(S_a * std::exp(S_a * x) * std::cos(S_a * y + S_d * z) - S_d * std::exp(S_a * z) * std::sin(S_a * x + S_d * y));
                case 1:
                    return -S_a * e *(S_a * std::exp(S_a * y) * std::sin(S_a * z + S_d * x) - S_a * std::exp(S_a * x) * std::sin(S_a * y + S_d * z));
                case 2:
                    return -S_a * e *(S_d * std::exp(S_a * z) * std::cos(S_a * x + S_d * y) + S_a * std::exp(S_a * y) * std::cos(S_a * z + S_d * x));
                default:
                    exit(1);
            }
        case 2:
            switch(i)
            {
                case 0:
                    return -S_a * e *(S_d * std::exp(S_a * x) * std::cos(S_a * y + S_d * z) +  S_a * std::exp(S_a * z) * std::cos(S_a * x + S_d * y));
                case 1:
                    return -S_a * e *(S_a * std::exp(S_a * y) * std::cos(S_a * z + S_d * x) - S_d * std::exp(S_a * x) * std::sin(S_a * y + S_d * z));
                case 2:
                    return -S_a * e *(S_a * std::exp(S_a * z) * std::sin(S_a * x + S_d * y) - S_a * std::exp(S_a * y) * std::sin(S_a * z + S_d * x));
                default:
                    exit(1);
            }
        default:
            exit(1);
    }
    return 1.;
}

double
RossEthierSteinmanDec::
grad_u_dt(const UInt& icoor, const double& t, const double& x, const double& y, const double& z, const ID& i)
{
    double e = -std::exp(-S_d * S_d * S_nu * t) * S_d * S_d * S_nu;
    switch(icoor)
    {
        case 0:    // u_x
            switch(i)
            {
                case 0:
                    return -S_a * e *(S_a * std::exp(S_a * x) * std::sin(S_a * y + S_d * z) - S_a * std::exp(S_a * z) * std::sin(S_a * x + S_d * y));
                case 1:
                    return -S_a * e *(S_d * std::exp(S_a * y) * std::cos(S_a * z + S_d * x) + S_a * std::exp(S_a * x) * std::cos(S_a * y + S_d * z));
                case 2:
                    return -S_a * e *(S_a * std::exp(S_a * z) * std::cos(S_a * x + S_d * y) - S_d * std::exp(S_a * y) * std::sin(S_a * z + S_d * x));
                default:
                    exit(1);
            }
        case 1:   // u_y
            switch(i)
            {
                case 0:
                    return -S_a * e *(S_a * std::exp(S_a * x) * std::cos(S_a * y + S_d * z) - S_d * std::exp(S_a * z) * std::sin(S_a * x + S_d * y));
                case 1:
                    return -S_a * e *(S_a * std::exp(S_a * y) * std::sin(S_a * z + S_d * x) - S_a * std::exp(S_a * x) * std::sin(S_a * y + S_d * z));
                case 2:
                    return -S_a * e *(S_d * std::exp(S_a * z) * std::cos(S_a * x + S_d * y) + S_a * std::exp(S_a * y) * std::cos(S_a * z + S_d * x));
                default:
                    exit(1);
            }
        case 2:
            switch(i)
            {
                case 0:
                    return -S_a * e *(S_d * std::exp(S_a * x) * std::cos(S_a * y + S_d * z) +  S_a * std::exp(S_a * z) * std::cos(S_a * x + S_d * y));
                case 1:
                    return -S_a * e *(S_a * std::exp(S_a * y) * std::cos(S_a * z + S_d * x) - S_d * std::exp(S_a * x) * std::sin(S_a * y + S_d * z));
                case 2:
                    return -S_a * e *(S_a * std::exp(S_a * z) * std::sin(S_a * x + S_d * y) - S_a * std::exp(S_a * y) * std::sin(S_a * z + S_d * x));
                default:
                    exit(1);
            }
        default:
            exit(1);
    }
    return 1.;
}

double
RossEthierSteinmanDec::f(const double& /* t */,
                                 const double& /* x */,
                                 const double& /* y */,
                                 const double& /* z */,
                                 const ID& /* i */)
{
    return 0.;
}

double
RossEthierSteinmanDec::
f_dt(const double& /* t */, const double& /* x */, const double& /* y */,
     const double& /* z */, const ID& i)
{
    return 0;
}

double
RossEthierSteinmanDec::xexact(const double& t,
                                      const double& x,
                                      const double& y,
                                      const double& z,
                                      const ID& i)
{
    switch(i)
    {
        case 0:
        case 1:
        case 2:
            return uexact(t, x, y, z, i);
        case 3:
            return pexact(t, x, y, z, 1);
        default:
            exit(1);
    }
    return 1.;
}



// Initial velocity
double
RossEthierSteinmanDec::x0(const double& t, const double& x, const double& y,
                                  const double& z, const ID& i)
{
    return xexact(t, x, y, z, i);
}

double
RossEthierSteinmanDec::fNeumann(const double& t,
                                const double& x,
                                const double& y,
                                const double& z,
                                const ID& i,
                                const VectorSmall<3>& n)
{
    double out = 0.;

    for(UInt k = 0; k < 3; k++) //mu grad_u n
    {
        out += S_mu * grad_u(k, t, x, y, z, i) * n[k];
    }

    if(S_flagStrain)
        for(UInt k = 0; k < 3; k++) //mu grad_u^T n
        {
            out += S_mu * grad_u(i, t, x, y, z, k) * n[k];
        }

    out -= pexact(t, x, y, z, i) * n[i]; //grad_p n
    return out;
}

double
RossEthierSteinmanDec::fNeumann_dt(const double& t,
                                   const double& x,
                                   const double& y,
                                   const double& z,
                                   const ID& i,
                                   const VectorSmall<3>& n)
{
    double out = 0.;

    for(UInt k = 0; k < 3; k++) //mu grad_u n
    {
        out += S_mu * grad_u_dt(k, t, x, y, z, i) * n[k];
    }

    if(S_flagStrain)
        for(UInt k = 0; k < 3; k++) //mu grad_u^T n
        {
            out += S_mu * grad_u_dt(i, t, x, y, z, k) * n[k];
        }

    out -= pexact_dt(t, x, y, z, i) * n[i]; //grad_p n
    return out;
}

void RossEthierSteinmanDec::setParamsFromGetPot(const GetPot& dataFile)
{
    S_a = dataFile("fluid/a", 1.) ;
    S_d = dataFile("fluid/d", 1.) ;
    S_mu = dataFile("fluid/viscosity", 1.);
    S_rho = dataFile("fluid/density", 1.);
    S_L = dataFile("fluid/L", 1.0);
    S_nu = S_mu / S_rho;
    S_flagStrain = dataFile("fluid/flag_strain", 0);
}

void
RossEthierSteinmanDec::
setA(const double& aValue)
{
    S_a = aValue;
}

void
RossEthierSteinmanDec::
setD(const double& dValue)
{
    S_d = dValue;
}

void
RossEthierSteinmanDec::
setViscosity(const double& mu)
{
    S_mu = mu;
    S_nu = S_mu / S_rho;
}

void
RossEthierSteinmanDec::
setDensity(const double& rho)
{
    S_rho = rho;
    S_nu = S_mu / S_rho;
}

void
RossEthierSteinmanDec::
setFlagStrain(const Int& flagValue)
{
    S_flagStrain = flagValue;
}

double RossEthierSteinmanDec::S_a;
double RossEthierSteinmanDec::S_d;
double RossEthierSteinmanDec::S_mu;
double RossEthierSteinmanDec::S_rho;
double RossEthierSteinmanDec::S_nu;
double RossEthierSteinmanDec::S_L;
int  RossEthierSteinmanDec::S_flagStrain;

} // namespace LifeV

#else

#include "RossEthierSteinmanDec.hpp"

namespace LifeV
{

double
RossEthierSteinmanDec::
uexact(const double& t, const double& x, const double& y,
       const double& z, const ID& i)
{
    switch(i)
    {
        case 0:
            return 0;
        case 1:
            return 0;
        case 2:
            return sin(x);//y * y * y;
        default:
            exit(1);
    }
    return 1.;
}

double
RossEthierSteinmanDec::
uexact_dt(const double& t, const double& x, const double& y,
          const double& z, const ID& i)
{
    double e = -std::exp(-S_d * S_d * S_nu * t) * S_d * S_d * S_nu;

    switch(i)
    {
        case 0:
            return 0;
        case 1:
            return 0;
        case 2:
            return 0;
        default:
            exit(1);
    }
    return 1.;
}

double
RossEthierSteinmanDec::
uderexact(const double& t, const double& x, const double& y,
          const double& z, const ID& i)
{

    if(i < 3)
    {
        return -S_d * S_d * S_nu * xexact(t, x, y, z, i);
    }
    else
    {
        return 0.;
    }
}

double
RossEthierSteinmanDec::
pexact(const double& t, const double& x, const double& y,
                                             const double& z,
                                             const ID& /* i */)
{
    return 0;
}

double
RossEthierSteinmanDec::
pexact_dt(const double& t, const double& x, const double& y,
                                             const double& z,
                                             const ID& /* i */)
{
    return 0;
}

double
RossEthierSteinmanDec::
grad_u(const UInt& icoor, const double& t, const double& x, const double& y, const double& z, const ID& i)
{
    double e = std::exp(-S_d * S_d * S_nu * t);
    switch(icoor)
    {
        case 0:    // u_x
            switch(i)
            {
                case 0:
                    return 0;
                case 1:
                    return 0;
                case 2:
                    return cos(x);
                default:
                    exit(1);
            }
        case 1:   // u_y
            switch(i)
            {
                case 0:
                    return 0;
                case 1:
                    return 0;
                case 2:
                    return 0;
                default:
                    exit(1);
            }
        case 2:
            switch(i)
            {
                case 0:
                    return 0;
                case 1:
                    return 0;
                case 2:
                    return 0;
                default:
                    exit(1);
            }
        default:
            exit(1);
    }
    return 1.;
}

double
RossEthierSteinmanDec::
grad_u_dt(const UInt& icoor, const double& t, const double& x, const double& y,
          const double& z, const ID& i)
{
    return 0.;
}

double
RossEthierSteinmanDec::
f(const double& /* t */, const double& x, const double& y,
  const double& /* z */, const ID& i)
{
    switch(i)
    {
        case 0:
            return 0;
        case 1:
            return 0;
        case 2:
            return sin(x);
        default:
            exit(1);
    }
    return 1.;
}

double
RossEthierSteinmanDec::
f_dt(const double& /* t */, const double& /* x */, const double& /* y */,
     const double& /* z */, const ID& i)
{
    switch(i)
    {
        case 0:
            return 0;
        case 1:
            return 0;
        case 2:
            return 0;
        default:
            exit(1);
    }
    return 1.;
}

double
RossEthierSteinmanDec::xexact(const double& t,
                                      const double& x,
                                      const double& y,
                                      const double& z,
                                      const ID& i)
{
    switch(i)
    {
        case 0:
        case 1:
        case 2:
            return uexact(t, x, y, z, i);
        case 3:
            return pexact(t, x, y, z, 1);
        default:
            exit(1);
    }
    return 1.;
}



// Initial velocity
double
RossEthierSteinmanDec::
x0(const double& t, const double& x, const double& y,
   const double& z, const ID& i)
{
    return xexact(t, x, y, z, i);
}



//we suppose that the problem geometry is the cube [-L,L]x[-L,L]x[-L,L].
double
RossEthierSteinmanDec::
fNeumann(const double& t, const double& x, const double& y,
         const double& z, const ID& i)
{

    double n[3] = {0., 0., 0.};
    double out = 0.;
    if (std::abs(x + S_L) < 1e-10)
    {
        n[0] = -1.;
    }
    else if(std::abs(x - S_L) < 1e-10)
    {
        n[0] =  1.;
    }
    else if(std::abs(y + S_L) < 1e-10)
    {
        n[1] = -1.;
    }
    else if(std::abs(x - S_L) < 1e-10)
    {
        n[1] =  1.;
    }
    else if(std::abs(z + S_L) < 1e-10)
    {
        n[2] = -1.;
    }
    else if(std::abs(z - S_L) < 1e-10)
    {
        n[2] =  1.;
    }
    else
    {
        std::cout << "strange point: x=" << x << " y=" << y << " z=" << z
                  << std::endl;
    }

    for(UInt k = 0; k < 3; k++) //mu grad_u n
    {
        out += S_mu * grad_u(k, t, x, y, z, i) * n[k];
    }

    if(S_flagStrain)
        for(UInt k = 0; k < 3; k++) //mu grad_u^T n
        {
            out += S_mu * grad_u(i, t, x, y, z, k) * n[k];
        }

    out -= pexact(t, x, y, z, i) * n[i]; //grad_p n
    // std::cout << "x = " << x << " y = " << y << " z = " << z << " i = " << i << std::endl << std::flush;
    // std::cout << out << std::endl << std::flush;
    return out;
}

double
RossEthierSteinmanDec::fNeumann_dt(const double& t,
                                           const double& x,
                                           const double& y,
                                           const double& z,
                                           const ID& i)
{

    double n[3] = {0., 0., 0.};
    double out = 0.;
    if (std::abs(x + S_L) < 1e-10)
    {
        n[0] = -1.;
    }
    else if(std::abs(x - S_L) < 1e-10)
    {
        n[0] =  1.;
    }
    else if(std::abs(y + S_L) < 1e-10)
    {
        n[1] = -1.;
    }
    else if(std::abs(x - S_L) < 1e-10)
    {
        n[1] =  1.;
    }
    else if(std::abs(z + S_L) < 1e-10)
    {
        n[2] = -1.;
    }
    else if(std::abs(z - S_L) < 1e-10)
    {
        n[2] =  1.;
    }
    else
    {
        std::cout << "strange point: x=" << x << " y=" << y << " z=" << z
                  << std::endl;
    }

    for(UInt k = 0; k < 3; k++) //mu grad_u n
    {
        out += S_mu * grad_u_dt(k, t, x, y, z, i) * n[k];
    }

    if(S_flagStrain)
        for(UInt k = 0; k < 3; k++) //mu grad_u^T n
        {
            out += S_mu * grad_u_dt(i, t, x, y, z, k) * n[k];
        }

    out -= pexact_dt(t, x, y, z, i) * n[i]; //grad_p n
    return out;
}

void RossEthierSteinmanDec::setParamsFromGetPot(const GetPot& dataFile)
{
    S_a = dataFile("fluid/a", 1.) ;
    S_d = dataFile("fluid/d", 1.) ;
    S_mu = dataFile("fluid/viscosity", 1.);
    S_rho = dataFile("fluid/density", 1.);
    S_L = dataFile("fluid/L", 1.0);
    S_nu = S_mu / S_rho;
    S_flagStrain = dataFile("fluid/flag_strain", 0);
}

void
RossEthierSteinmanDec::
setA(const double& aValue)
{
    S_a = aValue;
}

void
RossEthierSteinmanDec::
setD(const double& dValue)
{
    S_d = dValue;
}

void
RossEthierSteinmanDec::
setViscosity(const double& mu)
{
    S_mu = mu;
    S_nu = S_mu / S_rho;
}

void
RossEthierSteinmanDec::
setDensity(const double& rho)
{
    S_rho = rho;
    S_nu = S_mu / S_rho;
}

void
RossEthierSteinmanDec::
setFlagStrain(const Int& flagValue)
{
    S_flagStrain = flagValue;
}

double RossEthierSteinmanDec::S_a;
double RossEthierSteinmanDec::S_d;
double RossEthierSteinmanDec::S_mu;
double RossEthierSteinmanDec::S_rho;
double RossEthierSteinmanDec::S_nu;
double RossEthierSteinmanDec::S_L;
int  RossEthierSteinmanDec::S_flagStrain;

} // namespace LifeV
#endif
