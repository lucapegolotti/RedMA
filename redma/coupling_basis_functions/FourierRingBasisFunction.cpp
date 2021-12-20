#include "FourierRingBasisFunction.hpp"

namespace RedMA {

FourierRingBasisFunction::
FourierRingBasisFunction(const GeometricFace& face,
                         int nFrequenciesTheta,
                         const double mesh_size):
    BasisFunctionFunctor(face)
{
    M_nFrequenciesTheta = nFrequenciesTheta;

    if (M_nFrequenciesTheta < 0)
        M_nBasisFunctions = 0;
    else
    {
        M_nBasisFunctions = 2 * M_nFrequenciesTheta + 1;

        double pid2 = M_PI / 2;

        M_thetaFreq.push_back(0.0);
        M_thetaPhase.push_back(pid2);
        for (unsigned int i = 0; i < M_nFrequenciesTheta; i++)
        {
            M_thetaFreq.push_back((i + 1) * 2);
            M_thetaFreq.push_back((i + 1) * 2);
            M_thetaPhase.push_back(0.0);
            M_thetaPhase.push_back(pid2);
        }

        for (unsigned int i = 0; i < (2 * M_nFrequenciesTheta + 1); i++)
                M_auxIndicesTheta.push_back(i);

        // in this way I impose that, at R-h, the Gaussian is 1% of its maximum
        M_sigmaRadial = std::min(M_face.M_radius / 100.0,
                                 mesh_size / std::sqrt(2*std::log(100)));

        printlog(YELLOW, "Adding Fourier modes on the ring...");
    }

    M_type = "fourier_ring";
}

FourierRingBasisFunction::return_Type
FourierRingBasisFunction::
operator()(const Vector3D &pos) {
    double theta;
    double r;
    getThetaAndRadius(pos, theta, r);

    double R = M_face.M_radius;
    double scaleRadial = 1.0 / std::sqrt(2.0 * M_PI * std::pow(M_sigmaRadial, 2));

    unsigned int indexTheta = M_auxIndicesTheta[M_index];

    double returnVal =
            std::sin(M_thetaFreq[indexTheta] * theta + M_thetaPhase[indexTheta]) *
            scaleRadial * std::exp(-1.0 * std::pow(r - R, 2) / (2.0 * std::pow(M_sigmaRadial, 2)));

    return returnVal;
}

} // Namespace RedMA

