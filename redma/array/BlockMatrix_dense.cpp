#include "BlockMatrix.hpp"

namespace RedMA
{

// template <>
// template <>
// void
// BlockMatrix<DenseMatrix>::
// getRowProperty(std::shared_ptr<LifeV::MapEpetra>& outMap,
//                const unsigned int& rowIndex) const
// {
// }
//
// template <>
// template <>
// void
// BlockMatrix<DenseMatrix>::
// getColProperty(std::shared_ptr<LifeV::MapEpetra>& outMap,
//                const unsigned int& colIndex) const
// {
// }
//
// template <>
// template <>
// void
// BlockMatrix<DenseMatrix>::
// getRowsProperty(std::shared_ptr<LifeV::MapEpetra>& outMap) const
// {
// }
//
// template <>
// template <>
// void
// BlockMatrix<DenseMatrix>::
// getColsProperty(std::shared_ptr<LifeV::MapEpetra>& outMap) const
// {
// }

template <>
template <>
void
BlockMatrix<BlockMatrix<DenseMatrix>>::
convertVectorType(const DenseVector& inputVector,
                  BlockVector<BlockVector<DenseVector>>& outputVector) const
{
    if (!M_isFinalized)
        throw new Exception("Matrix must be finalized in order to convert vectors");

    outputVector.resize(M_nRows);

    unsigned int offset = 0;
    for (unsigned int i = 0; i < M_nRows; i++)
    {
        outputVector.block(i).resize(block(i,0).nRows());
        for (unsigned int ii = 0; ii < block(i,0).nRows(); ii++)
        {
            outputVector.block(i).block(ii).data().reset(new DENSEVECTOR(block(i,0).block(ii,0).getNumRows()));
            for (unsigned int iii = 0; iii < block(i,0).block(ii,0).getNumRows(); iii++)
            {
                (*outputVector.block(i).block(ii).data())(iii) = (*inputVector.data())(iii + offset);
            }
            offset += block(i,0).block(ii,0).getNumRows();
        }
    }
}

template <>
template <>
void
BlockMatrix<DenseMatrix>::
convertVectorType(const BlockVector<VectorEp>& inputVector,
                  Epetra_MultiVector& outputVector) const
{
}

template <>
void
BlockMatrix<DenseMatrix>::
finalize(BlockMatrix<BlockMatrix<DenseMatrix>>* father,
         unsigned int* myRow, unsigned int* myCol)
{
    std::vector<unsigned int> rows(M_nRows);
    std::vector<unsigned int> cols(M_nCols);

    if (father == nullptr)
    {
        for (unsigned int i = 0; i < M_nRows; i++)
        {
            for (unsigned int j = 0; j < M_nCols; j++)
            {
                if (block(i,j).getNumRows())
                    rows[i] = block(i,j).getNumRows();
                if (block(i,j).getNumCols())
                    cols[j] = block(i,j).getNumCols();
            }
        }

    }
    else
    {
        for (unsigned int i = 0; i < M_nRows; i++)
        {
            for (unsigned int curCol = 0; curCol < father->nCols(); curCol++)
            {
                unsigned int ncols = father->block(*myRow,curCol).nCols();
                for (unsigned int j = 0; j < ncols; j++)
                {
                    if (father->block(*myRow,curCol).block(i,j).getNumRows())
                    {
                        rows[i] = father->block(*myRow,curCol).block(i,j).getNumRows();
                        j = ncols;
                        curCol = father->nCols();
                    }
                }
            }
        }

        for (unsigned int j = 0; j < M_nCols; j++)
        {
            for (unsigned int curRow = 0; curRow < father->nRows(); curRow++)
            {
                unsigned int nrows = father->block(curRow,*myCol).nRows();
                for (unsigned int i = 0; i < nrows; i++)
                {
                    if (father->block(curRow,*myCol).block(i,j).getNumCols())
                    {
                        cols[j] = father->block(curRow,*myCol).block(i,j).getNumCols();
                        i = nrows;
                        curRow = father->nRows();
                    }
                }
            }
        }

    }

    for (unsigned int i = 0; i < M_nRows; i++)
    {
        for (unsigned int j = 0; j < M_nCols; j++)
        {
            if (block(i,j).isNull())
            {
                block(i,j).setNumRows(rows[i]);
                block(i,j).setNumCols(cols[j]);
            }
        }
    }

    M_isNull = true;
    for (unsigned int i = 0; i < M_nRows; i++)
    {
        for (unsigned int j = 0; j < M_nCols; j++)
        {
            M_isNull = M_isNull && block(i,j).isNull();
        }
    }

    M_isFinalized = true;
}

template <>
void
BlockMatrix<BlockMatrix<DenseMatrix>>::
finalize(BlockMatrix<BlockMatrix<BlockMatrix<DenseMatrix>>>* father,
         unsigned int* myRow, unsigned int* myCol)
{
    if (father != nullptr)
        throw new Exception("Error in finalize: BlockMatrix<BlockMatrix<BlockMatrix<DenseMatrix>>> not supported");


    std::vector<unsigned int> nrows(M_nRows);
    std::vector<unsigned int> ncols(M_nCols);

    for (unsigned int i = 0; i < M_nRows; i++)
    {
        for (unsigned int j = 0; j < M_nCols; j++)
        {
            if (block(i,j).nRows() > 0)
                nrows[i] = block(i,j).nRows();
            if (block(i,j).nCols() > 0)
                ncols[j] = block(i,j).nCols();
        }
    }

    for (unsigned int i = 0; i < M_nRows; i++)
    {
        for (unsigned int j = 0; j < M_nCols; j++)
        {
            if (block(i,j).nRows() == 0 && block(i,j).nCols() == 0)
                block(i,j).resize(nrows[i], ncols[j]);
            block(i,j).finalize(this, &i, &j);
        }
    }

    M_isNull = true;
    for (unsigned int i = 0; i < M_nRows; i++)
    {
        for (unsigned int j = 0; j < M_nCols; j++)
        {
            M_isNull = M_isNull && block(i,j).isNull();
        }
    }

    M_isFinalized = true;
}

template <>
void
BlockMatrix<DenseMatrix>::
printPattern() const
{
    for (unsigned int i = 0; i < M_nRows; i++)
    {
        for (unsigned int j = 0; j < M_nCols; j++)
        {
            auto curblock = block(i,j);
            if (curblock.data())
                std::cout << "(" << curblock.data()->M() <<
                             "," << curblock.data()->N() << ")";
            else
                std::cout << "o";
            std::cout << "\t\t";
        }
        std::cout << "\n";
    }
}

template <>
void
BlockMatrix<BlockMatrix<DenseMatrix>>::
printPattern() const
{
    for (unsigned int i = 0; i < M_nRows; i++)
    {
        for (unsigned int j = 0; j < M_nCols; j++)
        {
            auto curblock = block(i,j);
            if (!curblock.isNull())
                std::cout << "x";
            else
                std::cout << "o";
            std::cout << "\t";
        }
        std::cout << "\n";
    }
}

template <>
BlockMatrix<BlockMatrix<DenseMatrix>>
BlockMatrix<BlockMatrix<DenseMatrix>>::
getSubmatrix(const unsigned int& ibegin, const unsigned int& iend,
             const unsigned int& jbegin, const unsigned int& jend) const
{

    BlockMatrix<BlockMatrix<DenseMatrix>> retMatrix;

    if (!M_isFinalized)
        throw new Exception("Matrix must be finalized in order to apply getSubmatrix");

    unsigned int nrows = iend-ibegin+1;
    unsigned int ncols = jend-jbegin+1;
    retMatrix.resize(iend-ibegin+1, jend-jbegin+1);

    for (unsigned int i = ibegin; i <= iend; i++)
    {
        for (unsigned int j = jbegin; j <= jend; j++)
        {
            retMatrix.block(i-ibegin,j-jbegin).softCopy(block(i,j));
        }
    }

    retMatrix.finalize();

    return retMatrix;
}

template <>
DenseMatrix
BlockMatrix<DenseMatrix>::
collapse() const
{
    if (!M_isFinalized)
        throw new Exception("Matrix must be finalized in order to be collapsed");

    std::vector<unsigned int> nrows;
    std::vector<unsigned int> ncols;

    nrows.resize(M_nRows);
    ncols.resize(M_nCols);

    for (unsigned int i = 0; i < M_nRows; i++)
        nrows[i] = block(i,0).getNumRows();

    for (unsigned int j = 0; j < M_nCols; j++)
        ncols[j] = block(0,j).getNumCols();

    unsigned int totalrows = 0;
    for (auto row : nrows)
        totalrows += row;

    unsigned int totalcols = 0;
    for (auto col : ncols)
        totalcols += col;

    std::shared_ptr<DENSEMATRIX> inMatrix(new DENSEMATRIX(totalrows, totalcols));
    inMatrix->Scale(0.0);

    unsigned int offsetrow = 0;
    for (unsigned int i = 0; i < M_nRows; i++)
    {
        unsigned int offsetcol = 0;
        for (unsigned int j = 0; j < M_nCols; j++)
        {
            if (block(i,j).data())
            {
                for (unsigned ii = 0; ii < nrows[i]; ii++)
                {
                    for (unsigned jj = 0; jj < ncols[j]; jj++)
                    {
                        (*inMatrix)(ii + offsetrow, jj + offsetcol) = (*block(i,j).data())(ii,jj);
                    }
                }
            }
            offsetcol += ncols[j];
        }
        offsetrow += nrows[i];
    }

    DenseMatrix retMatrix;
    retMatrix.data() = inMatrix;

    return retMatrix;
}

// in principle we could call collapse recursively on the subblocks. However,
// I implement it from scratch in order to reduce the allocation of new matrices
// and to operate in place as much as possible
template <>
BlockMatrix<DenseMatrix>
BlockMatrix<BlockMatrix<DenseMatrix>>::
collapse() const
{
    if (!M_isFinalized)
        throw new Exception("Matrix must be finalized in order to be collapsed");

    std::vector<unsigned int> nrowsout;
    std::vector<unsigned int> ncolsout;
    std::vector<unsigned int> nrows;
    std::vector<unsigned int> ncols;

    for (unsigned int i = 0; i < M_nRows; i++)
    {
        unsigned int totrows = 0;
        for (unsigned int ii = 0; ii < block(i,0).nRows(); ii++)
        {
            totrows += block(i,0).block(ii,0).getNumRows();
            nrows.push_back(block(i,0).block(ii,0).getNumRows());
        }
        nrowsout.push_back(totrows);
    }

    for (unsigned int j = 0; j < M_nCols; j++)
    {
        unsigned int totcols = 0;
        for (unsigned int jj = 0; jj < block(0,j).nCols(); jj++)
        {
            totcols += block(0,j).block(0,jj).getNumCols();
            ncols.push_back(block(0,j).block(0,jj).getNumCols());
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
    for (unsigned int i = 0; i < M_nRows; i++)
    {
        unsigned int offsetcol;
        unsigned int offsetcolout = 0;
        for (unsigned int j = 0; j < M_nCols; j++)
        {
            auto curblock = block(i,j);
            unsigned int curRows = curblock.nRows();
            unsigned int curCols = curblock.nCols();

            unsigned int curCurRows;
            unsigned int curCurCols;

            offsetrow = offsetrowout;
            for (unsigned int ii = 0; ii < curRows; ii++)
            {
                offsetcol = offsetcolout;
                for (unsigned int jj = 0; jj < curCols; jj++)
                {
                    auto curmatrix = curblock.block(ii,jj);
                    curCurRows = curmatrix.getNumRows();
                    curCurCols = curmatrix.getNumCols();

                    if (curmatrix.data())
                    {
                        for (unsigned int iii = 0; iii < curCurRows; iii++)
                        {
                            for (unsigned int jjj = 0; jjj < curCurCols; jjj++)
                            {
                                (*inMatrix)(iii + offsetrow, jjj + offsetcol) =
                                              (*curmatrix.data())(iii,jjj);
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

    BlockMatrix<DenseMatrix> retMatrix(1,1);
    retMatrix.block(0,0).data() = inMatrix;

    return retMatrix;
}

}
