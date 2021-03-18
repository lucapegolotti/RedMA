#include "DummyBasisFunction.hpp"

namespace RedMA
{

DummyBasisFunction::
DummyBasisFunction(const GeometricFace& face, std::string type) :
  BasisFunctionFunctor(face)
{
    M_nBasisFunctions = 0;
    M_type = type;
}

DummyBasisFunction::return_Type
DummyBasisFunction::
operator()(const Vector3D& pos)
{
    // throw Exception("Dummy basis function should not be used!");
    return 0.0;
}

}  // namespace RedMA
