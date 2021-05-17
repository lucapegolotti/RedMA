#include "PreconditionerOperator.hpp"

namespace RedMA
{

PreconditionerOperator::
PreconditionerOperator()
{
}

void
PreconditionerOperator::
setPressureMass(const BM &mass)
{
    M_Mp.reset(new BlockMatrix());
    M_Mp->deepCopy(mass);
}

} // Namespace RedMA
