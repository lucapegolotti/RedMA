#include "ClothFunctionFunctor.hpp"

namespace RedMA {

ClothFunctionFunctor::
ClothFunctionFunctor(const LifeV::Vector3D &center, const double &radius,
                     const LifeV::Vector3D& normal, const LifeV::Vector3D& tangent,
                     const LifeV::Vector3D& shape_coefficients) :
M_center(center), M_radius(radius), M_normal(normal), M_tangent1(tangent), M_shapeCoeffs(shape_coefficients)
{
    M_normal.normalize();
    M_tangent1.normalize();

    M_tangent2 = M_normal.cross(M_tangent1);
    M_tangent2.normalize();

    // setting the eigenvector matrix
    Matrix3D eigenMatrix;

    eigenMatrix(0,0) = M_normal[0];
    eigenMatrix(0,1) = M_tangent1[0];
    eigenMatrix(0,2) = M_tangent2[0];
    eigenMatrix(1,0) = M_normal[1];
    eigenMatrix(1,1) = M_tangent1[1];
    eigenMatrix(1,2) = M_tangent2[1];
    eigenMatrix(2,0) = M_normal[2];
    eigenMatrix(2,1) = M_tangent1[2];
    eigenMatrix(2,2) = M_tangent2[2];

    // setting the eigenvalues matrix
    Matrix3D diagMatrix;

    diagMatrix(0,0) = shape_coefficients[0];
    diagMatrix(0,1) = 0.0;
    diagMatrix(0,2) = 0.0;
    diagMatrix(1,0) = 0.0;
    diagMatrix(1,1) = shape_coefficients[1];
    diagMatrix(1,2) = 0.0;
    diagMatrix(2,0) = 0.0;
    diagMatrix(2,1) = 0.0;
    diagMatrix(2,2) = shape_coefficients[2];

    // custom matrix norm
    M_normMatrix = eigenMatrix * diagMatrix * eigenMatrix.inverse();
}

ClothFunctionFunctor::Function
ClothFunctionFunctor::
function()
{
    return std::bind(&ClothFunctionFunctor::evaluateOperator, this,
                     std::placeholders::_1,
                     std::placeholders::_2,
                     std::placeholders::_3,
                     std::placeholders::_4,
                     std::placeholders::_5);
}

ClothFunctionFunctor::return_Type
ClothFunctionFunctor::
evaluateOperator(const double& t, const double& x, const double& y,
                 const double& z, unsigned int const& index)
{
    Vector3D pos(x,y,z);
    return this->operator()(pos);
}

ClothFunctionFunctor::return_Type
ClothFunctionFunctor::
operator()(const Vector3D& pos)
{
    Vector3D diff = pos - M_center;
    double dist = this->norm(diff);
    double bd_tol = 0.1;  // "tolerance" at the boundary

    if (dist <= (1.0 - bd_tol)*M_radius)
        return 1.0;
    else if ((dist > (1.0 - bd_tol)*M_radius) && (dist <= M_radius))
        return cos(M_PI*(dist - (1.0 - bd_tol)*M_radius) / (2*bd_tol*M_radius));
    else
        return 0.0;
}

ClothFunctionFunctor::return_Type
ClothFunctionFunctor::
norm(const LifeV::Vector3D &pos) const
{
    return pos.dot((M_normMatrix * pos));
}

}