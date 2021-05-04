#include "BasisFunctionFactory.hpp"

namespace RedMA
{

shp<BasisFunctionFunctor>
BasisFunctionFactory(const GetPot& datafile,
                     GeometricFace inlet,
                     bool isBoundary,
                     bool isRing)
{
    shp<BasisFunctionFunctor> basisFunction;
    std::string type = datafile("coupling/type", "chebyshev");

    if (isRing)
    {
        int frequenciesTheta = datafile("coupling/frequencies_theta_ring", -1);
        if (isBoundary)
            frequenciesTheta = datafile("coupling/frequencies_theta_ring_weak_dirichlet", -1);

        basisFunction.reset(new FourierRingBasisFunction(inlet,
                                                         frequenciesTheta));
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
            if (isBoundary)
                nMax = datafile("coupling/nfcts_weak_dirichlet", 1);
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
