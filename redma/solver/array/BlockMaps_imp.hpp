namespace RedMA
{

template <class InMatrixType>
BlockMaps<InMatrixType>::
BlockMaps(const BlockMatrix<InMatrixType>& matrix) :
  M_matrix(matrix)
{
    if (!M_matrix.isFinalized())
        throw new Exception("Matrix must be finalized before creating maps");
    generateMaps();
}

template <class InMatrixType>
SHP(MAPEPETRA)
BlockMaps<InMatrixType>::
getMonolithicRangeMapEpetra() const
{
    SHP(MAPEPETRA) retMap(new MAPEPETRA());

    for (auto map : M_rangeMapsEpetra)
        *retMap += *map;

    return retMap;
}

template <class InMatrixType>
SHP(MAPEPETRA)
BlockMaps<InMatrixType>::
getMonolithicDomainMapEpetra() const
{
    SHP(MAPEPETRA) retMap(new MAPEPETRA());

    for (auto map : M_domainMapsEpetra)
    {
        *retMap += *map;
    }

    return retMap;
}

}
