#include "BasisFunctionFunctor.hpp"

namespace RedMA
{

BasisFunctionFunctor::
BasisFunctionFunctor(const GeometricFace& face) :
  M_face(face)
{
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
    Vector3D& tangent1 = M_face.M_tangent1;
    Vector3D& tangent2 = M_face.M_tangent2;

    Vector3D diff = pos - center;

    x = diff.dot(tangent1);
    y = diff.dot(tangent2);
}

void
BasisFunctionFunctor::
getThetaAndRadius(const Vector3D& pos, double& theta, double& radius)
{
    Vector3D& center = M_face.M_center;
    Vector3D& normal = M_face.M_normal;
    Vector3D& tangent1 = M_face.M_tangent1;
    Vector3D& tangent2 = M_face.M_tangent2;

    Vector3D diff = pos - center;
    radius = diff.norm();

    if (radius < 1e-15)
        theta = 0;
    else
    {
        double ratio;
        if (diff.dot(tangent2) > 0)
        {
            ratio = diff.dot(tangent1) / radius;

            if (std::abs(ratio + 1) < 1e-15)
                theta = M_PI;
            else
                theta = std::acos(ratio);
        }
        else
        {
            ratio = -diff.dot(tangent1) / radius;

            if (std::abs(ratio + 1) < 1e-15)
                theta = M_PI;
            else
                theta = std::acos(ratio);
            theta += M_PI;
        }

    }
}

}  // namespace RedMA
