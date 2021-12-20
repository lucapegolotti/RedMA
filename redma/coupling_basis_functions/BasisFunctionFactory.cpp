#include "BasisFunctionFactory.hpp"

namespace RedMA
{

shp<BasisFunctionFunctor>
BasisFunctionFactory(const GetPot& datafile,
                     GeometricFace inlet,
                     bool isInlet, bool isOutlet,
                     bool isRing, const double mesh_size)
{
    shp<BasisFunctionFunctor> basisFunction;
    std::string type = datafile("coupling/type", "chebyshev");

    if (isRing)
    {
        int frequenciesTheta = datafile("coupling/frequencies_theta_ring", -1);
        if (isInlet)
            frequenciesTheta = datafile("coupling/frequencies_theta_ring_weak_dirichlet_in", -1);
        else if (isOutlet)
            frequenciesTheta = datafile("coupling/frequencies_theta_ring_weak_dirichlet_out", -1);

        basisFunction.reset(new FourierRingBasisFunction(inlet,
                                                           frequenciesTheta,
                                                           mesh_size));
    }

    else
    {
        if (!std::strcmp(type.c_str(), "fourier"))
        {
            unsigned int frequenciesTheta = datafile("coupling/frequencies_theta", 1);
            unsigned int frequenciesRadial = datafile("coupling/frequencies_radial", 1);

            basisFunction.reset(new FourierBasisFunction(inlet,
                                                           frequenciesTheta,
                                                           frequenciesRadial));
        }
        else if (!std::strcmp(type.c_str(), "chebyshev"))
        {
            unsigned int nMax = datafile("coupling/nMax", 5);
            if (isInlet)
                nMax = datafile("coupling/nfcts_weak_dirichlet_in", 5);
            else if (isOutlet)
                nMax = datafile("coupling/nfcts_weak_dirichlet_out", 0);
            basisFunction.reset(new ChebyshevBasisFunction(inlet, nMax));
        }
        else if (!std::strcmp(type.c_str(), "traces"))  // it actually defines null functions!
        {
            basisFunction.reset(new DummyBasisFunction(inlet, type));
        }
        else
        {
            std::string msg("Basis functions of type " + type + " not implemented!");
            throw Exception(msg);
        }
    }

    return basisFunction;
}

}  // namespace RedMA
