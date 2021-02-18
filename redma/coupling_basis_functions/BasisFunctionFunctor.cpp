#include "BasisFunctionFunctor.hpp"

namespace RedMA
{

BasisFunctionFunctor::
BasisFunctionFunctor(const GeometricFace& face) :
  M_face(face)
{
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

    M_eOrth = normal.cross(M_e);
    M_eOrth = M_eOrth / M_eOrth.norm();
}

void
BasisFunctionFunctor::
setIndex(const unsigned int& index)
{
    M_index = index;
}

unsigned int
BasisFunctionFunctor::
getNumBasisFunctions() const
{
    return M_nBasisFunctions;
}

BasisFunctionFunctor::Function
BasisFunctionFunctor::
function()
{
    return std::bind(&BasisFunctionFunctor::evaluateOperator, this,
                     std::placeholders::_1,
                     std::placeholders::_2,
                     std::placeholders::_3,
                     std::placeholders::_4,
                     std::placeholders::_5);
}

BasisFunctionFunctor::return_Type
BasisFunctionFunctor::
evaluateOperator(const double& t, const double& x, const double& y,
                 const double& z, unsigned int const& index)
{
    Vector3D pos(x,y,z);
    return this->operator()(pos);
}

void
BasisFunctionFunctor::
getLocalXAndY(const Vector3D& pos, double& x, double& y)
{
    Vector3D& center = M_face.M_center;

    Vector3D diff = pos - center;

    x = diff.dot(M_e);
    y = diff.dot(M_eOrth);
}

void
BasisFunctionFunctor::
getThetaAndRadius(const Vector3D& pos, double& theta, double& radius)
{
    Vector3D& center = M_face.M_center;
    Vector3D& normal = M_face.M_normal;

    Vector3D diff = pos - center;
    radius = diff.norm();

    if (radius < 1e-15)
        theta = 0;
    else
    {
        double ratio;
        if (diff.dot(M_eOrth) > 0)
        {
            ratio = diff.dot(M_e) / radius;

            if (std::abs(ratio + 1) < 1e-15)
                theta = M_PI;
            else
                theta = std::acos(ratio);
        }
        else
        {
            ratio = -diff.dot(M_e) / radius;

            if (std::abs(ratio + 1) < 1e-15)
                theta = M_PI;
            else
                theta = std::acos(ratio);
            theta += M_PI;
        }

    }
}

}  // namespace RedMA
