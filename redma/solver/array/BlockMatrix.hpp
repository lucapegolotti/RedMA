// Reduced Modeling of Arteries (RedMA)
// Copyright (C) 2019  Luca Pegolotti
//
// RedMA is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// RedMA is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with this program.  If not, see <http://www.gnu.org/licenses/>.

#ifndef BLOCKMATRIX_HPP
#define BLOCKMATRIX_HPP

#include <redma/utils/Exception.hpp>
#include <redma/solver/array/aMatrix.hpp>
#include <redma/solver/array/BlockVector.hpp>
#include <redma/solver/array/VectorEp.hpp>
#include <redma/solver/array/MatrixEp.hpp>

#include <boost/numeric/ublas/matrix.hpp>

#include <lifev/core/array/MapEpetra.hpp>


namespace RedMA
{

template <class InMatrixType>
class BlockMatrix : public aMatrix
{
    typedef boost::numeric::ublas::matrix<InMatrixType>       Grid;
public:

    typedef InMatrixType                                      InnerType;

    BlockMatrix();

    BlockMatrix(const unsigned int& nRows, const unsigned int& nCols);

    virtual BlockMatrix<InMatrixType> operator+(const BlockMatrix<InMatrixType>& other) const;

    virtual BlockMatrix<InMatrixType>& operator+=(const BlockMatrix<InMatrixType>& other);

    virtual BlockMatrix<InMatrixType>& operator*=(const double& coeff);

    virtual BlockMatrix<InMatrixType> operator*(const double& coeff) const;

    // hard copy!
    virtual BlockMatrix<InMatrixType>& operator=(const BlockMatrix<InMatrixType>& other);

    template <class InVectorType>
    BlockVector<InVectorType> operator*(const BlockVector<InVectorType>& vector) const;

    void resize(const unsigned int& nRows, const unsigned int& nCols);

    void hardCopy(const BlockMatrix<InMatrixType>& other);

    void softCopy(const BlockMatrix<InMatrixType>& other);

    InMatrixType& block(const unsigned int& iblock, const unsigned int& jblock);

    InMatrixType block(const unsigned int& iblock, const unsigned int& jblock) const;

    template <class InputVectorType, class OutputVectorType>
    void convertVectorType(const InputVectorType& inputVector,
                           OutputVectorType& output) const;

    // if InMatrix is epetra this returns global range map, otherwise number
    // of rows
    template <class OutputType>
    void getRowsProperty(OutputType& output) const;

    template <class OutputType>
    void getColsProperty(OutputType& output) const;

    template <class OutputType>
    void getRowProperty(OutputType& output, const unsigned int& indexrow) const;

    template <class OutputType>
    void getColProperty(OutputType& output, const unsigned int& indexcol) const;

    void collapseBlocks(InMatrixType& output);

    inline unsigned int nRows() const {return M_nRows;}

    inline unsigned int nCols() const {return M_nCols;}

    // this is to uniform all the dimensions in case of BlockMatrix<BlockMatrix<MatrixEp>>
    // (for now)
    void finalize();

protected:
    Grid            M_matrixGrid;
    unsigned int    M_nRows;
    unsigned int    M_nCols;
    // this is relevant only when block<block<matrix>> at the moment
    bool            M_isFinalized;
};

}

#include "BlockMatrix_imp.hpp"

#endif // BLOCKMATRIX_HPP
