#include <BasisFunctionFunctor.hpp>

namespace RedMA
{

BasisFunctionFunctor::
BasisFunctionFunctor(const GeometricFace& face) :
  M_face(face)
{
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
    using namespace std::placeholders;
    return std::bind(&BasisFunctionFunctor::evaluateOperator, this,
                     _1, _2, _3, _4, _5);
}

BasisFunctionFunctor::return_Type
BasisFunctionFunctor::
evaluateOperator(const double& t, const double& x, const double& y,
                 const double& z, unsigned int const& index)
{
    Vector3D pos(x,y,z);
    return this->operator()(pos);
}



}  // namespace RedMA
