#include <GlobalIdentityOperator.hpp>

#include <lifev/core/linear_algebra/IfpackPreconditioner.hpp>
#include <lifev/core/linear_algebra/MLPreconditioner.hpp>
#include <lifev/core/linear_algebra/TwoLevelPreconditioner.hpp>
#include <lifev/core/linear_algebra/AztecooOperatorAlgebra.hpp>
#include <lifev/core/linear_algebra/BelosOperatorAlgebra.hpp>
#include <lifev/core/linear_algebra/ApproximatedInvertibleRowMatrix.hpp>

#include <lifev/core/algorithm/PreconditionerML.hpp>

namespace LifeV
{
namespace Operators
{

GlobalIdentityOperator::
GlobalIdentityOperator() :
  M_label("GlobalIdentityOperator"),
  M_useTranspose(false)
{

}

GlobalIdentityOperator::~GlobalIdentityOperator()
{

}

void GlobalIdentityOperator::showMe(){}

void
GlobalIdentityOperator::
setUp(operatorPtrContainer_Type oper, const commPtr_Type& comm)
{
}


int GlobalIdentityOperator::ApplyInverse(const vector_Type& X, vector_Type& Y) const
{
    ASSERT_PRE(X.NumVectors() == Y.NumVectors(), "X and Y must have the same number of vectors");

    Y = X;

    return 0;
}
}
}
