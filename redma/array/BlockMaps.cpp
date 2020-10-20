#include "BlockMaps.hpp"

namespace RedMA
{

// template <>
// void
// BlockMaps<MatrixEp>::
// generateMaps()
// {
//     unsigned int nrows = M_matrix.nRows();
//     unsigned int ncols = M_matrix.nCols();
//
//     M_rangeMaps.resize(nrows);
//     M_domainMaps.resize(ncols);
//     M_rangeMapsEpetra.resize(nrows);
//     M_domainMapsEpetra.resize(ncols);
//     M_nInnerRows.push_back(nrows);
//     M_nInnerCols.push_back(ncols);
//
//     for (unsigned int i = 0; i < nrows; i++)
//     {
//         for (unsigned int j = 0; j < ncols; j++)
//         {
//             M_rangeMaps[i] = M_matrix.block(i,j).data()->rangeMap().map(LifeV::Unique);
//             M_domainMaps[j] = M_matrix.block(i,j).data()->domainMap().map(LifeV::Unique);
//             M_rangeMapsEpetra[i].reset(new MAPEPETRA(M_matrix.block(i,j).data()->rangeMap()));
//             M_domainMapsEpetra[j].reset(new MAPEPETRA(M_matrix.block(i,j).data()->domainMap()));
//         }
//     }
// }
//
// template <>
// void
// BlockMaps<BlockMatrix<MatrixEp>>::
// generateMaps()
// {
//     unsigned int nrows = M_matrix.nRows();
//     unsigned int ncols = M_matrix.nCols();
//
//     // first we deal with ranges
//     for (unsigned int i = 0; i < nrows; i++)
//     {
//         for (unsigned int j = 0; j < ncols; j++)
//         {
//             BlockMaps<MatrixEp> localMaps(M_matrix.block(i,j));
//             auto ranges = localMaps.getRangeMaps();
//
//             for (auto map : ranges)
//                 M_rangeMaps.push_back(map);
//
//             auto rangesEpetra = localMaps.getRangeMapsEpetra();
//             for (auto map : rangesEpetra)
//                 M_rangeMapsEpetra.push_back(map);
//
//             M_nInnerRows.push_back(localMaps.getInnerRows()[0]);
//             break;
//         }
//     }
//
//     for (unsigned int j = 0; j < ncols; j++)
//     {
//         for (unsigned int i = 0; i < nrows; i++)
//         {
//             BlockMaps<MatrixEp > localMaps(M_matrix.block(i,j));
//             auto ranges = localMaps.getDomainMaps();
//
//             for (auto map : ranges)
//                 M_domainMaps.push_back(map);
//
//             auto domainEpetra = localMaps.getDomainMapsEpetra();
//             for (auto map : domainEpetra)
//                 M_domainMapsEpetra.push_back(map);
//
//             M_nInnerCols.push_back(localMaps.getInnerCols()[0]);
//             break;
//         }
//     }
// }
//
// template <>
// BlockMatrix<MatrixEp>
// collapseBlocks(const BlockMatrix<BlockMatrix<MatrixEp>>& matrix,
//                const BlockMaps<BlockMatrix<MatrixEp>>& maps)
// {
//     if (!matrix.isFinalized())
//     {
//         throw new Exception("Matrix must be finalized in order to collapse blocks!");
//     }
//
//     auto rangeMaps = maps.getRangeMaps();
//     auto domainMaps = maps.getRangeMaps();
//
//     unsigned int nrows = rangeMaps.size();
//     unsigned int ncols = domainMaps.size();
//
//     std::vector<unsigned int> ninnerrows = maps.getInnerRows();
//     std::vector<unsigned int> ninnercols = maps.getInnerCols();
//
//     BlockMatrix<MatrixEp> retMatrix;
//     retMatrix.resize(nrows, ncols);
//
//     unsigned int offsetrow = 0;
//     for (unsigned int iouter = 0; iouter < matrix.nRows(); iouter++)
//     {
//         unsigned int offsetcol = 0;
//         for (unsigned int jouter = 0; jouter < matrix.nCols(); jouter++)
//         {
//             auto curmatrix = matrix.block(iouter,jouter);
//
//             for (unsigned int iinner = 0; iinner < ninnerrows[iouter]; iinner++)
//             {
//                 for (unsigned int jinner = 0; jinner < ninnercols[jouter]; jinner++)
//                 {
//                     retMatrix.block(offsetrow+iinner,offsetcol+jinner).
//                               softCopy(curmatrix.block(iinner,jinner));
//                 }
//             }
//             offsetcol += ninnercols[jouter];
//         }
//         offsetrow += ninnerrows[iouter];
//     }
//
//     retMatrix.finalize();
//
//     return retMatrix;
// }
//
// template <>
// MatrixEp
// collapseBlocks(const BlockMatrix<MatrixEp>& matrix,
//                const BlockMaps<MatrixEp>& maps)
// {
//     using namespace LifeV::MatrixEpetraStructuredUtility;
//
//     MatrixEp retMatrix;
//     auto rangeMaps = maps.getRangeMapsEpetra();
//     auto domainMaps = maps.getRangeMapsEpetra();
//
//     retMatrix.data().reset(new MATRIXEPETRA(*maps.getMonolithicRangeMapEpetra()));
//
//     std::vector<unsigned int> dimensionsRows;
//     std::vector<unsigned int> dimensionsCols;
//
//     for (auto map : rangeMaps)
//         dimensionsRows.push_back(map->mapSize());
//
//     for (auto map : domainMaps)
//         dimensionsCols.push_back(map->mapSize());
//
//     LifeV::MatrixBlockStructure structure;
//     structure.setBlockStructure(dimensionsRows,
//                                 dimensionsCols);
//
//     for (unsigned int i = 0; i < matrix.nRows(); i++)
//     {
//         for (unsigned int j = 0; j < matrix.nCols(); j++)
//         {
//             SHP(LifeV::MatrixEpetraStructuredView<double>) globalView;
//             globalView = createBlockView(retMatrix.data(), structure, i, j);
//
//             LifeV::MatrixBlockStructure blockStructure;
//             std::vector<unsigned int> rows(1), cols(1);
//             rows[0] = dimensionsRows[i];
//             cols[0] = dimensionsCols[j];
//             blockStructure.setBlockStructure(rows, cols);
//
//             SHP(LifeV::MatrixEpetraStructuredView<double>) blockLocalView;
//             blockLocalView = createBlockView(matrix.block(i,j).data(),
//                                              blockStructure, 0, 0);
//
//             copyBlock(blockLocalView, globalView);
//         }
//     }
//
//     retMatrix.data()->globalAssemble(maps.getMonolithicDomainMapEpetra(),
//                                      maps.getMonolithicRangeMapEpetra());
//
//     return retMatrix;
// }
//
// SHP(VECTOREPETRA)
// getEpetraVector(const BlockVector<BlockVector<VectorEp>>& vector,
//                 const BlockMaps<BlockMatrix<MatrixEp>>& maps)
// {
//
//     SHP(VECTOREPETRA) retVec(new VECTOREPETRA(maps.getMonolithicRangeMapEpetra(),
//                                               LifeV::Unique));
//     auto rangeMaps = maps.getRangeMapsEpetra();
//     auto rows = maps.getInnerRows();
//
//     unsigned int count = 0;
//     unsigned int offset = 0;
//     // this part is tricky because we don't know a priori if all the
//     // blocks in the input vector are filled
//     for (unsigned int iouter = 0; iouter < vector.nRows(); iouter++)
//     {
//         auto curblock = vector.block(iouter);
//         unsigned int currows = vector.block(iouter).nRows();
//
//         if (curblock.nRows() == 0)
//         {
//             for (unsigned int iinner = 0; iinner < rows[iouter]; iinner++)
//             {
//                 offset += rangeMaps[count]->mapSize();
//                 count++;
//             }
//         }
//         else
//         {
//             for (unsigned int iinner = 0; iinner < curblock.nRows(); iinner++)
//             {
//                 if (curblock.block(iinner).data())
//                 {
//                     retVec->subset(*curblock.block(iinner).data(),
//                                    *rangeMaps[count], 0, offset);
//                 }
//
//                 offset += rangeMaps[count]->mapSize();
//                 count++;
//             }
//
//         }
//
//     }
//
//     return retVec;
// }
//
// BlockVector<BlockVector<VectorEp>>
// getBlockVector(const SHP(VECTOREPETRA)& vector,
//                const BlockMaps<BlockMatrix<MatrixEp>>& maps)
// {
//     BlockVector<BlockVector<VectorEp>> retVec;
//     std::vector<unsigned int> rows = maps.getInnerRows();
//     auto rangeMaps = maps.getRangeMapsEpetra();
//
//     retVec.resize(rows.size());
//
//     unsigned int count = 0;
//     unsigned int offset = 0;
//     for (unsigned int iouter = 0; iouter < rows.size(); iouter++)
//     {
//         retVec.block(iouter).resize(rows[iouter]);
//         for (unsigned int iinner = 0; iinner < rows[iouter]; iinner++)
//         {
//             SHP(VECTOREPETRA)& cvHandler = retVec.block(iouter).block(iinner).data();
//             cvHandler.reset(new VECTOREPETRA(rangeMaps[count], LifeV::Unique));
//             cvHandler->subset(*vector, *rangeMaps[count], offset, 0);
//             offset += rangeMaps[count]->mapSize();
//             count++;
//         }
//     }
//
//     return retVec;
// }

void
BlockMaps::
updateCollapsedMatrix(SHP(BlockMatrix) matrix)
{
    M_collapsedMatrix = collapseBlocks(matrix, M_dimensionsRowsBlock, M_dimensionsColsBlock);
}

void
BlockMaps::
createFromBlockMatrix(SHP(BlockMatrix) matrix)
{
    if (matrix->level() >= 2)
        updateCollapsedMatrix(matrix);
    else if (matrix->level() == 1)
        M_collapsedMatrix = matrix;
    matrix = M_collapsedMatrix;

    M_rangeEpetraMaps.resize(matrix->nRows(),nullptr);
    M_domainEpetraMaps.resize(matrix->nCols(),nullptr);
    M_rangeMaps.resize(matrix->nRows(),nullptr);
    M_domainMaps.resize(matrix->nCols(),nullptr);
    // collect maps and build monolithic
    for (unsigned int i = 0; i < matrix->nRows(); i++)
    {
        for (unsigned int j = 0; j < matrix->nCols(); j++)
        {
            if (matrix->block(i,j)->type() == SPARSE)
            {
                SHP(MATRIXEPETRA) curMatrix = std::static_pointer_cast<MATRIXEPETRA>(matrix->block(i,j)->data());
                if (curMatrix)
                {
                    M_rangeEpetraMaps[i] = curMatrix->rangeMap().map(LifeV::Unique);
                    M_domainEpetraMaps[j] = curMatrix->domainMap().map(LifeV::Unique);
                    M_rangeMaps[i].reset(new MAPEPETRA(curMatrix->rangeMap()));
                    M_domainMaps[j].reset(new MAPEPETRA(curMatrix->domainMap()));
                }
            }
        }
    }

    M_monolithicRangeMap.reset(new MAPEPETRA());
    for (unsigned int i = 0; i < matrix->nRows(); i++)
    {
        if (M_rangeMaps[i] == nullptr)
            throw new Exception("Error in createFromBlockMatrix!");

        *M_monolithicRangeMap += *M_rangeMaps[i];
        M_dimensionsRows.push_back(M_rangeMaps[i]->mapSize());
    }

    M_monolithicDomainMap.reset(new MAPEPETRA());
    for (unsigned int j = 0; j < matrix->nCols(); j++)
    {
        if (M_domainMaps[j] == nullptr)
            throw new Exception("Error in createFromBlockMatrix!");

        *M_monolithicDomainMap += *M_domainMaps[j];
        M_dimensionsCols.push_back(M_domainMaps[j]->mapSize());
    }
}


SHP(BlockMatrix)
collapseBlocks(SHP(BlockMatrix) matrix,
               std::vector<unsigned int>& dimensionsRowsBlock,
               std::vector<unsigned int>& dimensionsColsBlock)
{
    if (matrix->level() == 1)
        return matrix;

    if (matrix->level() == 2)
    {
        dimensionsRowsBlock.resize(matrix->nRows(),0);
        dimensionsColsBlock.resize(matrix->nCols(),0);

        for (unsigned int i = 0; i < matrix->nRows(); i++)
        {
            for (unsigned int j = 0; j < matrix->nCols(); j++)
            {
                if (!matrix->block(i,j)->isZero())
                {
                    if (matrix->block(i,j)->type() != BLOCK)
                        throw new Exception("Non BLOCK block in level-2 MatrixBlock!");
                    else
                    {
                        unsigned int nlocalrows = matrix->block(i,j)->nRows();
                        if (dimensionsRowsBlock[i] != 0 && dimensionsRowsBlock[i] != nlocalrows)
                            throw new Exception("Non consistent dimensions in BlockMatrix");
                        else
                            dimensionsRowsBlock[i] = nlocalrows;

                        unsigned int nlocalcols = matrix->block(i,j)->nCols();
                        if (dimensionsColsBlock[j] != 0 && dimensionsColsBlock[j] != nlocalcols)
                            throw new Exception("Non consistent dimensions in BlockMatrix");
                        else
                            dimensionsColsBlock[j] = nlocalcols;
                    }
                }
            }
        }

        std::vector<unsigned int> cumulativeRows(matrix->nRows(),0);

        unsigned int count = 0;
        unsigned int totalrows = 0;
        for (auto row : dimensionsRowsBlock)
        {
            totalrows += row;
            cumulativeRows[count] = totalrows;
            count++;
        }

        std::vector<unsigned int> cumulativeCols(matrix->nCols(),0);

        count = 0;
        unsigned int totalcols = 0;
        for (auto col : dimensionsColsBlock)
        {
            totalcols += col;
            cumulativeCols[count] = totalcols;
            count++;
        }


        SHP(BlockMatrix) retMat(new BlockMatrix(totalrows,totalcols));
        for (unsigned int i = 0; i < totalrows; i++)
        {
            for (unsigned int j = 0; j < totalcols; j++)
            {
                SHP(aMatrix) newBlock;

                unsigned int indexrow = 0;
                while (cumulativeRows[indexrow] < i + 1)
                    indexrow++;

                unsigned int indexcol = 0;
                while (cumulativeCols[indexcol] < j + 1)
                    indexcol++;

                if (matrix->block(indexrow,indexcol)->isZero())
                {
                    newBlock.reset(new DenseMatrix());
                }
                else
                {
                    unsigned int localrow = (indexrow == 0) ? i : i - cumulativeRows[indexrow-1];
                    unsigned int localcol = (indexcol == 0) ? j : j - cumulativeCols[indexcol-1];
                    newBlock = matrix->block(indexrow,indexcol)->block(localrow, localcol);
                }
                retMat->setBlock(i,j,newBlock);
            }
        }

        return retMat;
    }

    if (matrix->level() >= 3)
        throw new Exception ("Case level >= 3 not yet implemented!");


    return nullptr;
}

SHP(SparseMatrix)
blockMatrixToSparseMatrix(SHP(BlockMatrix) matrix)
{
    using namespace LifeV::MatrixEpetraStructuredUtility;

    BlockMaps maps(matrix);

    matrix = maps.M_collapsedMatrix;

    SHP(MATRIXEPETRA) ptrMatrix(new MATRIXEPETRA(*maps.M_monolithicRangeMap));
    //
    // std::vector<unsigned int> dimensionsRows;
    // std::vector<unsigned int> dimensionsCols;
    //
    // for (auto map : rangeMaps)
    //     dimensionsRows.push_back(map->mapSize());
    //
    // for (auto map : domainMaps)
    //     dimensionsCols.push_back(map->mapSize());
    //
    LifeV::MatrixBlockStructure structure;
    structure.setBlockStructure(maps.M_dimensionsRows,
                                maps.M_dimensionsCols);

    for (unsigned int i = 0; i < matrix->nRows(); i++)
    {
        for (unsigned int j = 0; j < matrix->nCols(); j++)
        {
            SHP(LifeV::MatrixEpetraStructuredView<double>) globalView;
            globalView = createBlockView(ptrMatrix, structure, i, j);

            LifeV::MatrixBlockStructure blockStructure;
            std::vector<unsigned int> rows(1), cols(1);
            rows[0] = maps.M_dimensionsRows[i];
            cols[0] = maps.M_dimensionsCols[j];
            blockStructure.setBlockStructure(rows, cols);

            SHP(LifeV::MatrixEpetraStructuredView<double>) blockLocalView;
            if (!matrix->block(i,j)->isZero())
            {
                blockLocalView = createBlockView(std::static_pointer_cast<MATRIXEPETRA>(matrix->block(i,j)->data()),
                                                 blockStructure, 0, 0);
                copyBlock(blockLocalView, globalView);
            }
        }
    }
    ptrMatrix->globalAssemble(maps.M_monolithicDomainMap, maps.M_monolithicRangeMap);

    SHP(SparseMatrix) retMatrix(new SparseMatrix());
    retMatrix->setMatrix(ptrMatrix);

    return retMatrix;
}

// void
// BlockDimension::
// close()
// {
//     // for (auto dim : M_dimensions)
//     // {
//     //     if (dim == nullptr)
//     //         throw new Exception("Error: one or more dimensions in BlockDimension is not set");
//     //
//     //     if (dim->primaryType() == BLOCK)
//     //         std::static_pointer_cast<BlockDimension>(dim)->close();
//     // }
//     //
//     // DimensionType first = M_dimensions[0]->secondaryType();
//     //
//     // for (auto dim : M_dimensions)
//     // {
//     //     if (dim->secondaryType() != first)
//     //     {
//     //         setSecondaryType(MIXED);
//     //         break;
//     //     }
//     // }
//     //
//     // if (secondaryType() == BLOCK)
//     // {
//     //     M_
//     // }
// }

SHP(aVector)
getBlockVector(const SHP(VECTOREPETRA)& vector, const BlockMaps& maps)
{
    std::vector<unsigned int> rows = maps.M_dimensionsRowsBlock;
    auto rangeMaps = maps.M_rangeMaps;

    SHP(BlockVector) retVec(new BlockVector(rows.size()));

    unsigned int count = 0;
    unsigned int offset = 0;
    for (unsigned int iouter = 0; iouter < rows.size(); iouter++)
    {
        SHP(BlockVector) curBlock(new BlockVector(rows[iouter]));
        for (unsigned int iinner = 0; iinner < rows[iouter]; iinner++)
        {
            SHP(VECTOREPETRA) innerBlockPtr(new VECTOREPETRA(rangeMaps[count], LifeV::Unique));
            innerBlockPtr->subset(*vector, *rangeMaps[count], offset, 0);
            offset += rangeMaps[count]->mapSize();
            count++;
            SHP(DistributedVector) distrVector(new DistributedVector());
            distrVector->setVector(innerBlockPtr);
            curBlock->setBlock(iinner, distrVector);
        }
        retVec->setBlock(iouter, curBlock);
    }

    return retVec;
}

SHP(VECTOREPETRA)
getEpetraVector(const SHP(aVector)& vector, const BlockMaps& maps)
{
    SHP(VECTOREPETRA) retVec(new VECTOREPETRA(maps.M_monolithicRangeMap,
                                              LifeV::Unique));
    auto rangeMaps = maps.M_rangeMaps;
    auto rows = maps.M_dimensionsRowsBlock;

    unsigned int count = 0;
    unsigned int offset = 0;
    // this part is tricky because we don't know a priori if all the
    // blocks in the input vector are filled
    for (unsigned int iouter = 0; iouter < vector->nRows(); iouter++)
    {
        auto curblock = vector->block(iouter);
        unsigned int currows = vector->block(iouter)->nRows();

        if (curblock->nRows() == 0 || curblock->isZero())
        {
            for (unsigned int iinner = 0; iinner < rows[iouter]; iinner++)
            {
                offset += rangeMaps[count]->mapSize();
                count++;
            }
        }
        else
        {
            for (unsigned int iinner = 0; iinner < curblock->nRows(); iinner++)
            {
                if (curblock->block(iinner)->data())
                {
                    SHP(VECTOREPETRA) localVector;
                    if (curblock->block(iinner)->type() == DENSE)
                    {
                        localVector = std::static_pointer_cast<DenseVector>(curblock->block(iinner))->toVectorEpetraPtr(rangeMaps[count]->commPtr());
                    }
                    else
                    {
                        localVector = std::static_pointer_cast<VECTOREPETRA>(curblock->block(iinner)->data());
                    }
                    retVec->subset(*localVector,
                                   *rangeMaps[count], 0, offset);
                }

                offset += rangeMaps[count]->mapSize();
                count++;
            }
        }
    }
    return retVec;
}

// SHP(VECTOREPETRA)
// getEpetraVector(const SHP(aVector)& vector, const BlockMaps& maps)
// {
//     std::cout << "getEpetraVector" << std::endl << std::flush;
//     SHP(VECTOREPETRA) retVec(new VECTOREPETRA(maps.M_monolithicRangeMap,
//                                               LifeV::Unique));
//     auto rangeMaps = maps.M_rangeMaps;
//     auto rows = maps.M_dimensionsRowsBlock;
//
//     std::cout << "1" << std::endl << std::flush;
//     unsigned int count = 0;
//     unsigned int offset = 0;
//     // this part is tricky because we don't know a priori if all the
//     // blocks in the input vector are filled
//     for (unsigned int iouter = 0; iouter < vector->nRows(); iouter++)
//     {
//         std::cout << "2" << std::endl << std::flush;
//         auto curblock = vector->block(iouter);
//         unsigned int currows = vector->block(iouter)->nRows();
//
//         if (curblock->nRows() == 0 || curblock->isZero())
//         {
//             for (unsigned int iinner = 0; iinner < rows[iouter]; iinner++)
//             {
//                 offset += rangeMaps[count]->mapSize();
//                 count++;
//             }
//         }
//         else
//         {
//             std::cout << "3" << std::endl << std::flush;
//             for (unsigned int iinner = 0; iinner < curblock->nRows(); iinner++)
//             {
//                 std::cout << "4" << std::endl << std::flush;
//                 if (curblock->block(iinner)->data())
//                 {
//                     std::cout << "4_" << std::endl << std::flush;
//                     SHP(VECTOREPETRA) localVector;
//                     if (curblock->block(iinner)->type() == DENSE)
//                     {
//                         std::cout << "4_1" << std::endl << std::flush;
//                         localVector = std::static_pointer_cast<DenseVector>(curblock->block(iinner))->toVectorEpetraPtr(rangeMaps[count]->commPtr());
//                     }
//                     else
//                     {
//                         std::cout << "4_2" << std::endl << std::flush;
//                         localVector = std::static_pointer_cast<VECTOREPETRA>(curblock->block(iinner)->data());
//                     }
//                     std::cout << "4_3" << std::endl << std::flush;
//                     std::cout << localVector->mapPtr()->mapSize() << std::endl << std::flush;
//                     std::cout << rangeMaps[count]->mapSize() << std::endl << std::flush;
//                     retVec->subset(*localVector,
//                                    *rangeMaps[count], 0, offset);
//                     std::cout << "4_4" << std::endl << std::flush;
//                 }
//                 std::cout << "5" << std::endl << std::flush;
//
//                 offset += rangeMaps[count]->mapSize();
//                 count++;
//                 std::cout << "6" << std::endl << std::flush;
//             }
//         }
//         std::cout << "7" << std::endl << std::flush;
//
//     }
//     return retVec;
// }

// SHP(DenseMatrix)
// blockMatrixToDenseMatrix(SHP(BlockMatrix) matrix)
// {
//
//
// }

}
