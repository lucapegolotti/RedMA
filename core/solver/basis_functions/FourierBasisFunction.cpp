#include <FourierBasisFunction.hpp>

namespace RedMA
{

FourierBasisFunction::
FourierBasisFunction(const GeometricFace& face,
                     unsigned int nFrequenciesTheta,
                     unsigned int nFrequenciesRadial) :
  BasisFunctionFunctor(face)
{
    M_nFrequenciesTheta = nFrequenciesTheta;
    M_nFrequenciesRadial = nFrequenciesRadial;
    // in the radial direction we only consider cos in order to ensure bfs
    // with continuous derivatives in the center
    M_nBasisFunctions = (2 * nFrequenciesTheta + 1) * (nFrequenciesRadial + 1);

    double pid2 = M_PI / 2;

    double radius = face.M_radius;

    M_thetaFreq.push_back(0.0);
    M_thetaPhase.push_back(pid2);
    for (unsigned int i = 0; i < nFrequenciesTheta; i++)
    {
        M_thetaFreq.push_back((i + 1) * 2);
        M_thetaFreq.push_back((i + 1) * 2);
        M_thetaPhase.push_back(0.0);
        M_thetaPhase.push_back(pid2);
    }

    M_radialPhase.push_back(0.0);
    M_radialFreq.push_back(0.0);
    for (unsigned int i = 0; i < nFrequenciesRadial; i++)
    {
        M_radialFreq.push_back((i + 1) * radius * 2);
        M_radialPhase.push_back(0.0);
    }

    for (unsigned int i = 0; i < (2 * nFrequenciesTheta + 1); i++)
    {
        for (unsigned int j = 0; j < nFrequenciesRadial + 1; j++)
        {
            M_auxIndicesTheta.push_back(i);
            M_auxIndicesRadial.push_back(j);
        }
    }
    M_type = "fourier";
}

FourierBasisFunction::return_Type
FourierBasisFunction::
operator()(const Vector3D& pos)
{
    double theta;
    double r;

    getThetaAndRadius(pos, theta, r);

    unsigned int indexTheta = M_auxIndicesTheta[M_index];
    unsigned int indexRadial = M_auxIndicesRadial[M_index];

    // we arbitrarily scale the basis functions in order to enduce an ordering
    // in the POD
    double scale = 1; // 1./ (std::pow(2, indexTheta + indexRadial));
    double returnVal =
           std::sin(M_thetaFreq[indexTheta] * theta + M_thetaPhase[indexTheta]) *
           std::cos(M_radialFreq[indexRadial] * r + M_radialPhase[indexRadial]) *
           scale;
    return returnVal;
}

}  // namespace RedMA
