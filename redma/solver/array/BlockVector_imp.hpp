#include "BlockVector.hpp"

namespace RedMA
{

template <class InVectorType>
BlockVector<InVectorType>::
BlockVector()
{
}

template <class InVectorType>
BlockVector<InVectorType>::
BlockVector(const unsigned int& nRows) :
  BlockMatrix<InVectorType>(nRows, 1)
{
}

template <class InVectorType>
InVectorType&
BlockVector<InVectorType>::
block(const unsigned int& iblock)
{
    return BlockMatrix<InVectorType>::block(iblock, 0);
}

template <class InVectorType>
InVectorType
BlockVector<InVectorType>::
block(const unsigned int& iblock) const
{
    return BlockMatrix<InVectorType>::block(iblock,0);
}

}
