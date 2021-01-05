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
#include <redma/array/aMatrix.hpp>
#include <redma/array/BlockVector.hpp>
#include <redma/array/DistributedVector.hpp>
#include <redma/array/DenseVector.hpp>
#include <redma/array/SparseMatrix.hpp>
#include <redma/array/DenseMatrix.hpp>
#include <redma/array/Double.hpp>
#include <redma/array/Wrap.hpp>

#include <boost/numeric/ublas/matrix.hpp>

#include <lifev/core/array/MapEpetra.hpp>


namespace RedMA
{

class BlockMatrix : public aMatrix
{

    typedef boost::numeric::ublas::matrix<shp<aMatrix>>       Grid;

public:

    BlockMatrix();

    virtual ~BlockMatrix() {};

    BlockMatrix(const unsigned int& nRows, const unsigned int& nCols);

    virtual void add(shp<aMatrix> other) override;

    virtual void multiplyByScalar(const double& coeff) override;

    virtual shp<aMatrix> multiplyByMatrix(shp<aMatrix> other) override;

    virtual shp<aMatrix> transpose() const override;

    virtual shp<aVector> multiplyByVector(shp<aVector> vector) override;

    virtual void shallowCopy(shp<aDataWrapper> other) override;

    virtual void deepCopy(shp<aDataWrapper> other) override;

    virtual bool isZero() const override;

    virtual void dump(std::string filename) const override;

    virtual BlockMatrix* clone() const override;

    void resize(const unsigned int& nRows, const unsigned int& nCols);

    shp<aMatrix> block(const unsigned int& iblock, const unsigned int& jblock) const;

    void setBlock(const unsigned int& iblock, const unsigned int& jblock,
                  shp<aMatrix> matrix);

    shp<BlockMatrix> getSubmatrix(const unsigned int& ibegin, const unsigned int& iend,
                                  const unsigned int& jbegin, const unsigned int& jend) const;

    unsigned int level();

    bool globalTypeIs(Datatype type);

    shp<BlockMatrix> convertInnerTo(Datatype type, shp<Epetra_Comm> comm = nullptr);

    void printPattern() const;

    void findComm();

    shp<Epetra_Comm> commPtr() {findComm(); return M_comm;}

    virtual Datatype type() const override {return BLOCK;}

    virtual shp<void> data() const override {return nullptr;};

    virtual void setData(shp<void> data) override {throw new Exception("setData undefined for BlockMatrix");};

protected:

    shp<Epetra_Comm>              M_comm;
    Grid                          M_matrixGrid;
};

shp<SparseMatrix> blockMatrixToSparseMatrix(shp<BlockMatrix> matrix);

}

#endif // BLOCKMATRIX_HPP
