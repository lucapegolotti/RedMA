#include "BasisFunctionFunctor.hpp"

namespace RedMA
{

BasisFunctionFunctor::
BasisFunctionFunctor(const GeometricFace& face) :
  M_face(face)
{
    // tangents are not considered in non-affine transformations --> handle tangents here
    if (std::abs(face.M_normal.dot(face.M_tangent1)) > 1e-8)
    {
        printlog(YELLOW, "[BasisFunctionFunctor] Invalid tangent vector at face with flag " +
                 std::to_string(face.M_flag) + ". Computing new tangent vectors. This is likely to yield "
                 "inaccurate results with elliptic faces!", true);

        if (std::abs(std::abs(face.M_normal[0]) - 1.0) > 1e-12)
        {
            M_e[0] = 1.0; M_e[1] = 0.0; M_e[2] = 0.0;
        }
        else
        {
            M_e(0) = 0.0; M_e[1] = 1.0; M_e[2] = 0.0;
        }
        // project the vector onto the face and orthonormalize
        M_e = M_e - M_e.dot(face.M_normal) * face.M_normal;
        M_e = M_e / M_e.norm();

        M_eOrth = face.M_normal.cross(M_e);
        M_eOrth = M_eOrth / M_eOrth.norm();
    }

    else
    {
        M_e = face.M_tangent1;
        M_eOrth = face.M_tangent2;
    }

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
