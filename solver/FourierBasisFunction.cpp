#include <FourierBasisFunction.hpp>

namespace RedMA
{

FourierBasisFunction::
FourierBasisFunction(const GeometricFace& face, unsigned int nFrequencies) :
  BasisFunctionFunctor(face)
{
    M_nFrequencies = nFrequencies;

    // in the radial direction we only consider cos in order to ensure bfs
    // with continuous derivatives in the center
    M_nBasisFunctions = (2 * nFrequencies + 1) * (nFrequencies + 1);

    double pid2 = M_PI / 2;

    double radius = face.M_radius;

    M_thetaFreq.push_back(0.0);
    M_thetaPhase.push_back(pid2);
    for (unsigned int i = 0; i < nFrequencies; i++)
    {
        M_thetaFreq.push_back(i);
        M_thetaFreq.push_back(i);
        M_thetaPhase.push_back(0.0);
        M_thetaPhase.push_back(pid2);
    }

    M_radialPhase.push_back(0.0);
    for (unsigned int i = 0; i < nFrequencies; i++)
    {
        M_radialFreq.push_back(i * radius);
        M_radialPhase.push_back(0.0);
    }

    for (unsigned int i = 0; i < (2 * nFrequencies + 1); i++)
    {
        for (unsigned int j = 0; j < nFrequencies + 1; j++)
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
        M_e[0] = 1.0; M_e[0] = 0.0; M_e[0] = 0.0;
    }
    else
    {
        M_e[0] = 0.0; M_e[0] = 1.0; M_e[0] = 0.0;
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

    if (crossProduct.dot(normal) < 0)
        theta = theta + M_PI;

    return std::sin(M_thetaFreq[M_index] * theta + M_thetaPhase[M_index]) *
           std::sin(M_radialFreq[M_index] * r + M_radialPhase[M_index]);
}

}  // namespace RedMA
