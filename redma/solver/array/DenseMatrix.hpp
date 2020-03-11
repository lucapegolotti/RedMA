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

#ifndef DENSEMATRIX_HPP
#define DENSEMATRIX_HPP


#include <Epetra_ConfigDefs.h>
#ifdef EPETRA_MPI
#include <mpi.h>
#include <Epetra_MpiComm.h>
#else
#include <Epetra_SerialComm.h>
#endif

#include <redma/utils/PrintLog.hpp>
#include <redma/solver/array/aMatrix.hpp>
#include <redma/solver/array/DenseVector.hpp>

#include <Epetra_SerialDenseMatrix.h>

#include <fstream>

#define DENSEMATRIX         Epetra_SerialDenseMatrix

namespace RedMA
{

class DenseMatrix : public aMatrix
{
public:
    DenseMatrix();

    // DenseMatrix(const std::vector<int>& columnVectors);

    DenseMatrix operator+(const DenseMatrix& other);

    DenseMatrix transpose() const;

    DenseMatrix& operator+=(const DenseMatrix& other);

    DenseMatrix& operator-=(const DenseMatrix& other);

    DenseMatrix& operator*=(const double& coeff);

    DenseVector operator*(const DenseVector& vector);

    DenseMatrix operator*(const DenseMatrix& other);

    void hardCopy(const DenseMatrix& other);

    void softCopy(const DenseMatrix& other);

    std::shared_ptr<DENSEMATRIX>& data();

    std::shared_ptr<DENSEMATRIX> data() const;

    void dump(std::string filename) const;

    // this method calls the corresponding method on the internal pointer.
    // If it is null, it returns the M_nRows
    unsigned int getNumRows() const;

    // this method calls the corresponding method on the internal pointer.
    // If it is null, it returns the M_nCols
    unsigned int getNumCols() const;

    void setNumRows(unsigned int numRows);

    void setNumCols(unsigned int numCols);

    inline bool isNull() const
    {
        if (!M_matrix)
            return true;
        return M_matrix->NormInf() < 1e-15;
    };

private:
    std::shared_ptr<DENSEMATRIX>  M_matrix;
    // differently from sparse matrices, we introduce num rows and cols. This is
    // because we don't want to allocate a dense matrix a zeros. So if in a block
    // system a matrix is zero, we store no pointer and we set the number of cols
    // rows
    unsigned int                  M_nCols;
    unsigned int                  M_nRows;
};

}

#endif // DENSEMATRIX_HPP
