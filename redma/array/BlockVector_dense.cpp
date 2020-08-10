#include "BlockVector.hpp"

namespace RedMA
{

template<>
DenseVector
BlockVector<DenseVector>::
collapse() const
{
    std::vector<unsigned int> nrows(M_nRows);

    unsigned int totalrows = 0;
    for (unsigned int i = 0; i < M_nRows; i++)
    {
        nrows[i] = block(i).getNumRows();
        totalrows += nrows[i];
    }

    std::shared_ptr<DENSEVECTOR> vector(new DENSEVECTOR(totalrows));

    unsigned int offset = 0;
    for (unsigned int i = 0; i < M_nRows; i++)
    {
        for (unsigned int ii = 0; ii < nrows[i]; ii++)
            (*vector)(ii + offset) = (*block(i).data())(ii);

        offset += nrows[i];
    }

    DenseVector retVector;
    retVector.data() = vector;

    return retVector;
}

template<>
BlockVector<DenseVector>
BlockVector<BlockVector<DenseVector>>::
collapse() const
{
    std::vector<unsigned int> nrowsout;
    std::vector<unsigned int> nrows;

    unsigned int totalrows = 0;
    for (unsigned int i = 0; i < M_nRows; i++)
    {
        unsigned int totrows = 0;
        for (unsigned int ii = 0; ii < block(i).nRows(); ii++)
        {
            totrows += block(i).block(ii).getNumRows();
            nrows.push_back(block(i).block(ii).getNumRows());
        }
        nrowsout.push_back(totrows);
        totalrows += totrows;
    }
    std::shared_ptr<DENSEVECTOR> vector(new DENSEVECTOR(totalrows));

    unsigned int offset = 0;
    for (unsigned int i = 0; i < M_nRows; i++)
    {
        for (unsigned int ii = 0; ii < block(i).nRows(); ii++)
        {
            for (unsigned int iii = 0; iii < block(i).block(ii).getNumRows(); iii++)
            {
                (*vector)(iii + offset) = (*block(i).block(ii).data())(iii);
            }
            offset += block(i).block(ii).getNumRows();
        }
    }
    BlockVector<DenseVector> retVector(1);
    retVector.block(0).data() = vector;
    return retVector;
}

}
