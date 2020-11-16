#include "BlockMaps.hpp"

namespace RedMA
{

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

}
