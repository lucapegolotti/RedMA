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
#include <redma/array/aMatrix.hpp>
#include <redma/array/DenseVector.hpp>

#include <Epetra_SerialDenseMatrix.h>

#include <fstream>

namespace RedMA
{

class DenseMatrix : public aMatrix
{
public:
    DenseMatrix();

    // DenseMatrix(const std::vector<int>& columnVectors);

    virtual void add(shp<aMatrix> other) override;

    virtual void multiplyByScalar(const double& coeff) override;

    virtual shp<aMatrix> multiplyByMatrix(shp<aMatrix> other) override;

    virtual shp<aMatrix> transpose() const override;

    virtual shp<aVector> multiplyByVector(shp<aVector> vector) override;

    virtual void shallowCopy(shp<aDataWrapper> other) override;

    virtual void deepCopy(shp<aDataWrapper> other) override;

    virtual bool isZero() const override;

    virtual shp<void> data() const override;

    virtual void setData(shp<void> data) override;

    virtual void dump(std::string filename) const override;

    virtual DenseMatrix* clone() const override;

    virtual Datatype type() const override {return DENSE;}

    void setMatrix(shp<DENSEMATRIX> matrix);

private:
    shp<DENSEMATRIX>        M_matrix;
};

}

#endif // DENSEMATRIX_HPP
