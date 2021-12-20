#include "BlockMaps.hpp"

namespace RedMA
{

void
BlockMaps::
updateCollapsedMatrix(shp<BlockMatrix> matrix)
{
    M_collapsedMatrix = collapseBlocks(matrix, M_dimensionsRowsBlock, M_dimensionsColsBlock);
}

void
BlockMaps::
createFromBlockMatrix(shp<BlockMatrix> matrix)
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
                shp<MATRIXEPETRA> curMatrix = spcast<MATRIXEPETRA>(matrix->block(i,j)->data());
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


shp<BlockMatrix>
collapseBlocks(shp<BlockMatrix> matrix,
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


        shp<BlockMatrix> retMat(new BlockMatrix(totalrows,totalcols));
        for (unsigned int i = 0; i < totalrows; i++)
        {
            for (unsigned int j = 0; j < totalcols; j++)
            {
                shp<aMatrix> newBlock;

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
                    newBlock = convert<BlockMatrix>(matrix->block(indexrow,indexcol))->block(localrow, localcol);
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

shp<SparseMatrix>
blockMatrixToSparseMatrix(shp<BlockMatrix> matrix)
{
    using namespace LifeV::MatrixEpetraStructuredUtility;

    BlockMaps maps(matrix);

    matrix = maps.M_collapsedMatrix;

    shp<MATRIXEPETRA> ptrMatrix(new MATRIXEPETRA(*maps.M_monolithicRangeMap));
    LifeV::MatrixBlockStructure structure;
    structure.setBlockStructure(maps.M_dimensionsRows,
                                maps.M_dimensionsCols);

    for (unsigned int i = 0; i < matrix->nRows(); i++)
    {
        for (unsigned int j = 0; j < matrix->nCols(); j++)
        {
            shp<LifeV::MatrixEpetraStructuredView<double>> globalView;
            globalView = createBlockView(ptrMatrix, structure, i, j);

            LifeV::MatrixBlockStructure blockStructure;
            std::vector<unsigned int> rows(1), cols(1);
            rows[0] = maps.M_dimensionsRows[i];
            cols[0] = maps.M_dimensionsCols[j];
            blockStructure.setBlockStructure(rows, cols);

            shp<LifeV::MatrixEpetraStructuredView<double>> blockLocalView;
            if (!matrix->block(i,j)->isZero())
            {
                blockLocalView = createBlockView(spcast<MATRIXEPETRA>(matrix->block(i,j)->data()),
                                                 blockStructure, 0, 0);
                copyBlock(blockLocalView, globalView);
            }
        }
    }
    ptrMatrix->globalAssemble(maps.M_monolithicDomainMap, maps.M_monolithicRangeMap);

    shp<SparseMatrix> retMatrix(new SparseMatrix());
    retMatrix->setMatrix(ptrMatrix);

    return retMatrix;
}

shp<DenseMatrix>
blockMatrixToDenseMatrix(shp<BlockMatrix> matrix)
{
    shp<DenseMatrix> retMatrix(new DenseMatrix());
    if (matrix->level() == 1)
    {
        std::vector<unsigned int> nrows;
        std::vector<unsigned int> ncols;

        unsigned int nRows = matrix->nRows();
        unsigned int nCols = matrix->nCols();

        nrows.resize(matrix->nRows());
        ncols.resize(matrix->nCols());

        for (unsigned int i = 0; i < nRows; i++)
            nrows[i] = matrix->block(i,0)->nRows();

        for (unsigned int j = 0; j < nCols; j++)
            ncols[j] = matrix->block(0,j)->nCols();

        unsigned int totalrows = 0;
        for (auto row : nrows)
            totalrows += row;

        unsigned int totalcols = 0;
        for (auto col : ncols)
            totalcols += col;

        std::shared_ptr<DENSEMATRIX> inMatrix(new DENSEMATRIX(totalrows, totalcols));
        inMatrix->Scale(0.0);

        unsigned int offsetrow = 0;
        for (unsigned int i = 0; i < nRows; i++)
        {
            unsigned int offsetcol = 0;
            for (unsigned int j = 0; j < nCols; j++)
            {
                if (matrix->block(i,j)->data())
                {
                    for (unsigned ii = 0; ii < nrows[i]; ii++)
                    {
                        for (unsigned jj = 0; jj < ncols[j]; jj++)
                        {
                            auto innerMatrix = spcast<DENSEMATRIX>(matrix->block(i,j)->data());
                            (*inMatrix)(ii + offsetrow, jj + offsetcol) = (*innerMatrix)(ii,jj);
                        }
                    }
                }
                offsetcol += ncols[j];
            }
            offsetrow += nrows[i];
        }

        retMatrix->setMatrix(inMatrix);
    }
    else if (matrix->level() == 2)
    {
        std::vector<unsigned int> nrowsout;
        std::vector<unsigned int> ncolsout;
        std::vector<unsigned int> nrows;
        std::vector<unsigned int> ncols;

        for (unsigned int i = 0; i < matrix->nRows(); i++)
        {
            unsigned int totrows = 0;
            unsigned int j = 0;
            while (matrix->block(i,j)->nRows() == 0)
                j++;
            for (unsigned int ii = 0; ii < matrix->block(i,j)->nRows(); ii++)
            {
                unsigned int jj = 0;
                while (matrix->block(i,j)->block(ii,jj)->nRows() == 0)
                    jj++;
                totrows += matrix->block(i,j)->block(ii,jj)->nRows();
                nrows.push_back(matrix->block(i,j)->block(ii,jj)->nRows());
            }
            nrowsout.push_back(totrows);
        }

        for (unsigned int j = 0; j < matrix->nCols(); j++)
        {
            unsigned int totcols = 0;
            unsigned int i = 0;
            while (matrix->block(i,j)->nCols() == 0)
                i++;
            for (unsigned int jj = 0; jj < matrix->block(i,j)->nCols(); jj++)
            {
                unsigned int ii = 0;
                while (matrix->block(i,j)->block(ii,jj)->nCols() == 0)
                    ii++;
                totcols += matrix->block(i,j)->block(ii,jj)->nCols();
                ncols.push_back(matrix->block(i,j)->block(ii,jj)->nCols());
            }
            ncolsout.push_back(totcols);
        }

        unsigned int totalrows = 0;
        for (auto row : nrows)
            totalrows += row;

        unsigned int totalcols = 0;
        for (auto col : ncols)
            totalcols += col;

        std::shared_ptr<DENSEMATRIX> inMatrix(new DENSEMATRIX(totalrows, totalcols));
        inMatrix->Scale(0.0);

        unsigned int offsetrow;
        unsigned int offsetrowout = 0;
        for (unsigned int i = 0; i < matrix->nRows(); i++)
        {
            unsigned int offsetcol;
            unsigned int offsetcolout = 0;
            for (unsigned int j = 0; j < matrix->nCols(); j++)
            {
                auto curblock = matrix->block(i,j);
                unsigned int curRows = curblock->nRows();
                unsigned int curCols = curblock->nCols();

                unsigned int curCurRows;
                unsigned int curCurCols;

                offsetrow = offsetrowout;
                for (unsigned int ii = 0; ii < curRows; ii++)
                {
                    offsetcol = offsetcolout;
                    for (unsigned int jj = 0; jj < curCols; jj++)
                    {
                        auto curmatrix = curblock->block(ii,jj);
                        curCurRows = curmatrix->nRows();
                        curCurCols = curmatrix->nCols();

                        if (curmatrix->data())
                        {
                            for (unsigned int iii = 0; iii < curCurRows; iii++)
                            {
                                for (unsigned int jjj = 0; jjj < curCurCols; jjj++)
                                {
                                    (*inMatrix)(iii + offsetrow, jjj + offsetcol) =
                                                  (*spcast<DENSEMATRIX>(curmatrix->data()))(iii,jjj);
                                }
                            }
                        }
                        offsetcol += curCurCols;
                    }
                    offsetrow += curCurRows;
                }
                offsetcolout += ncolsout[j];
            }
            offsetrowout += nrowsout[i];
        }

        retMatrix->setMatrix(inMatrix);
    }
    else
        throw new Exception("blockMatrixToDenseMatrix is only implemented for level <= 2");

    return retMatrix;
}

shp<aVector>
getBlockVector(const shp<VECTOREPETRA>& vector, const BlockMaps& maps)
{
    std::vector<unsigned int> rows = maps.M_dimensionsRowsBlock;
    auto rangeMaps = maps.M_rangeMaps;

    shp<BlockVector> retVec(new BlockVector(rows.size()));

    unsigned int count = 0;
    unsigned int offset = 0;
    for (unsigned int iouter = 0; iouter < rows.size(); iouter++)
    {
        shp<BlockVector> curBlock(new BlockVector(rows[iouter]));
        for (unsigned int iinner = 0; iinner < rows[iouter]; iinner++)
        {
            shp<VECTOREPETRA> innerBlockPtr(new VECTOREPETRA(rangeMaps[count], LifeV::Unique));
            innerBlockPtr->subset(*vector, *rangeMaps[count], offset, 0);
            offset += rangeMaps[count]->mapSize();
            count++;
            shp<DistributedVector> distrVector(new DistributedVector());
            distrVector->setVector(innerBlockPtr);
            curBlock->setBlock(iinner, distrVector);
        }
        retVec->setBlock(iouter, curBlock);
    }

    return retVec;
}

shp<VECTOREPETRA>
getEpetraVector(const shp<aVector>& vector, const BlockMaps& maps)
{
    shp<VECTOREPETRA> retVec(new VECTOREPETRA(maps.M_monolithicRangeMap,
                                              LifeV::Unique));
    auto rangeMaps = maps.M_rangeMaps;
    auto rows = maps.M_dimensionsRowsBlock;

    shp<BlockVector> blckVector = convert<BlockVector>(vector);

    unsigned int count = 0;
    unsigned int offset = 0;
    // this part is tricky because we don't know a priori if all the
    // blocks in the input vector are filled
    for (unsigned int iouter = 0; iouter < vector->nRows(); iouter++)
    {
        if (blckVector->block(iouter)->nRows() == 0 || blckVector->block(iouter)->isZero())
        {
            for (unsigned int iinner = 0; iinner < rows[iouter]; iinner++)
            {
                offset += rangeMaps[count]->mapSize();
                count++;
            }
        }
        else
        {
            auto curblock = convert<BlockVector>(blckVector->block(iouter));

            for (unsigned int iinner = 0; iinner < curblock->nRows(); iinner++)
            {
                if (curblock->block(iinner)->data())
                {
                    shp<VECTOREPETRA> localVector;
                    if (curblock->block(iinner)->type() == DENSE)
                    {
                        localVector = spcast<DenseVector>(curblock->block(iinner))->toVectorEpetraPtr(rangeMaps[count]->commPtr());
                    }
                    else
                    {
                        localVector = spcast<VECTOREPETRA>(curblock->block(iinner)->data());
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

shp<DenseVector>
blockVectorToDenseVector(shp<BlockVector> vector)
{
    shp<DenseVector> retVector(new DenseVector());

    if (vector->level() == 1)
    {
        std::vector<unsigned int> nrows(vector->nRows());

        unsigned int totalrows = 0;
        for (unsigned int i = 0; i < vector->nRows(); i++)
        {
            nrows[i] = vector->block(i)->nRows();
            totalrows += nrows[i];
        }

        std::shared_ptr<DENSEVECTOR> inVector(new DENSEVECTOR(totalrows));

        unsigned int offset = 0;
        for (unsigned int i = 0; i < vector->nRows(); i++)
        {
            for (unsigned int ii = 0; ii < nrows[i]; ii++)
                (*inVector)(ii + offset) = vector->block(i)->operator()(ii);

            offset += nrows[i];
        }
        retVector->setVector(inVector);
    }
    else if (vector->level() == 2)
    {
        std::vector<unsigned int> nrowsout;
        std::vector<unsigned int> nrows;

        unsigned int totalrows = 0;
        for (unsigned int i = 0; i < vector->nRows(); i++)
        {
            unsigned int totrows = 0;
            for (unsigned int ii = 0; ii < vector->block(i)->nRows(); ii++)
            {
                totrows += vector->block(i)->block(ii)->nRows();
                nrows.push_back(vector->block(i)->block(ii)->nRows());
            }
            nrowsout.push_back(totrows);
            totalrows += totrows;
        }
        std::shared_ptr<DENSEVECTOR> inVector(new DENSEVECTOR(totalrows));

        unsigned int offset = 0;
        for (unsigned int i = 0; i < vector->nRows(); i++)
        {
            for (unsigned int ii = 0; ii < vector->block(i)->nRows(); ii++)
            {
                for (unsigned int iii = 0; iii < vector->block(i)->block(ii)->nRows(); iii++)
                {
                    (*inVector)(iii + offset) = vector->block(i)->block(ii)->operator()(iii);
                }
                offset += vector->block(i)->block(ii)->nRows();
            }
        }
        retVector->setVector(inVector);
    }
    return retVector;
}

}
