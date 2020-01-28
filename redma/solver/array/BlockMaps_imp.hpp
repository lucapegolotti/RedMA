namespace RedMA
{

template <class InMatrixType>
BlockMaps<InMatrixType>::
BlockMaps(const BlockMatrix<InMatrixType>& matrix) :
  M_matrix(matrix)
{
    generateMaps();
}

}
