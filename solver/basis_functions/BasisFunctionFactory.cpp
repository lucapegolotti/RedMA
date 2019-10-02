#include <BasisFunctionFactory.hpp>

namespace RedMA
{

std::shared_ptr<BasisFunctionFunctor>
BasisFunctionFactory(const GetPot& datafile,
                     GeometricFace inlet)
{
    std::shared_ptr<BasisFunctionFunctor> basisFunction;
    std::string type = datafile("coupling/type", "chebyshev");

    if (!std::strcmp(type.c_str(), "fourier"))
    {
        unsigned int frequenciesTheta = datafile("coupling/frequencies_theta", 1);
        unsigned int frequenciesRadial = datafile("coupling/frequencies_radial", 1);

        basisFunction.reset(new FourierBasisFunction(inlet,
                                                     frequenciesTheta,
                                                     frequenciesRadial));
    }
    else if (!std::strcmp(type.c_str(), "zernike"))
    {
        unsigned int nMax = datafile("coupling/nMax", 5);
        basisFunction.reset(new ZernikeBasisFunction(inlet, nMax));
    }
    else if (!std::strcmp(type.c_str(), "chebyshev"))
    {
        unsigned int nMax = datafile("coupling/nMax", 5);
        basisFunction.reset(new ChebyshevBasisFunction(inlet, nMax));
    }
    else if (!std::strcmp(type.c_str(), "traces"))
    {
        basisFunction.reset(new DummyBasisFunction(inlet, type));
    }
    else
    {
        std::string msg("Basis functions of type " + type + " not implemented!");
        throw Exception(msg);
    }

    return basisFunction;
}

}  // namespace RedMA
