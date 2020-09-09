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

#ifndef SPARSEMATRIX_HPP
#define SPARSEMATRIX_HPP

#include <Epetra_ConfigDefs.h>
#ifdef EPETRA_MPI
#include <mpi.h>
#include <Epetra_MpiComm.h>
#else
#include <Epetra_SerialComm.h>
#endif

#include <redma/utils/PrintLog.hpp>
#include <redma/array/aMatrix.hpp>
#include <redma/array/DistributedVector.hpp>
#include <redma/array/DenseMatrix.hpp>

#include <lifev/core/array/MatrixEpetra.hpp>
#include <lifev/core/array/MapEpetra.hpp>

#define MATRIXEPETRA        LifeV::MatrixEpetra<double>

namespace RedMA
{

class SparseMatrix : public aMatrix
{
public:
    SparseMatrix();

    ~SparseMatrix() {};

    SparseMatrix(const SparseMatrix& other);

    SparseMatrix(const std::vector<DistributedVector>& columnVectors);

    SparseMatrix(const std::vector<std::shared_ptr<VECTOREPETRA>>& columnVectors);

    virtual void add(std::shared_ptr<aMatrix> other) override;

    virtual void multiplyByScalar(const double& coeff) override;

    virtual std::shared_ptr<aMatrix> multiplyByMatrix(std::shared_ptr<aMatrix> other) override;

    virtual std::shared_ptr<aMatrix> transpose() const override;

    virtual std::shared_ptr<aVector> multiplyByVector(std::shared_ptr<aVector> vector) override;

    virtual void softCopy(std::shared_ptr<aMatrix> other) override;

    virtual void hardCopy(std::shared_ptr<aMatrix> other) override;

    virtual bool isZero() const override;

    virtual SparseMatrix* clone() const override;

    void setMatrix(std::shared_ptr<MATRIXEPETRA> matrix);

    std::shared_ptr<void> data() const override;

    DenseMatrix toDenseMatrix() const;

    virtual void dump(std::string namefile) const override;

private:
    std::shared_ptr<MATRIXEPETRA>           M_matrix;
};

}

#endif // SPARSEMATRIX_HPP
