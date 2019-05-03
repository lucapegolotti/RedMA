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


}  // namespace RedMA
