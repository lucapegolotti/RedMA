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
        M_thetaFreq.push_back(i + 1);
        M_thetaFreq.push_back(i + 1);
        M_thetaPhase.push_back(0.0);
        M_thetaPhase.push_back(pid2);
    }

    M_radialPhase.push_back(0.0);
    M_radialFreq.push_back(0.0);
    for (unsigned int i = 0; i < nFrequenciesRadial; i++)
    {
        M_radialFreq.push_back((i + 1) * radius);
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

    Vector3D& normal = M_face.M_normal;

    // arbitrary vector to measure the angle
    // we check just the first component of the normal because we trust that
    // it is unitary (hence if normal[1] == 1 => normal = (1,0,0))
    if (std::abs(std::abs(normal[0]) - 1.0) > 1e-12)
    {
        M_e[0] = 1.0; M_e[1] = 0.0; M_e[2] = 0.0;
    }
    else
    {
        M_e[0] = 0.0; M_e[1] = 1.0; M_e[2] = 0.0;
    }
    // project the vector onto the face and orthonormalize
    M_e = M_e - M_e.dot(normal) * normal;
    M_e = M_e / M_e.norm();
}

FourierBasisFunction::return_Type
FourierBasisFunction::
operator()(const Vector3D& pos)
{
    Vector3D& center = M_face.M_center;
    Vector3D& normal = M_face.M_normal;

    Vector3D diff = pos - center;
    double r = diff.norm();

    double theta = std::acos(diff.dot(M_e) / (diff.norm()));
    Vector3D crossProduct = diff.cross(M_e);

    unsigned int indexTheta = M_auxIndicesTheta[M_index];
    unsigned int indexRadial = M_auxIndicesRadial[M_index];

    if (crossProduct.dot(normal) < 0)
        theta = theta + M_PI;
    double returnVal =
           std::sin(M_thetaFreq[indexTheta] * theta + M_thetaPhase[indexTheta]) *
           std::cos(M_radialFreq[indexRadial] * r + M_radialPhase[indexRadial]);
    return returnVal;
}

}  // namespace RedMA
