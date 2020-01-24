#include "LinearOperatorEp.hpp"

namespace RedMA
{

LinearOperatorEp::
LinearOperatorEp()
{
}

void
LinearOperatorEp::
setup(BM matrix, EPETRACOMM comm)
{
    M_matrix.softCopy(matrix);
    M_comm = comm;
}

int
LinearOperatorEp::
Apply(const super::vector_Type& X, super::vector_Type& Y) const
{
    BlockVector<VectorEp> xblock;
    M_matrix.convertVectorType(X, xblock);

    BlockVector<VectorEp> yblock = M_matrix * xblock;
    M_matrix.convertVectorType(yblock, Y);
}

}
