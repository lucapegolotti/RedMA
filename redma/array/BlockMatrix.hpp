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
#include <redma/array/SparseMatrix.hpp>
#include <redma/array/DenseMatrix.hpp>
#include <redma/array/Double.hpp>

#include <boost/numeric/ublas/matrix.hpp>

#include <lifev/core/array/MapEpetra.hpp>


namespace RedMA
{

class BlockMatrix : public aMatrix
{
    typedef boost::numeric::ublas::matrix<std::shared_ptr<aMatrix>>       Grid;
public:

    virtual ~BlockMatrix() {};

    BlockMatrix(const BlockMatrix& other);

    BlockMatrix(const unsigned int& nRows, const unsigned int& nCols);

    virtual void add(std::shared_ptr<aMatrix> other) override;

    virtual void multiplyByScalar(const double& coeff) override;

    virtual std::shared_ptr<aMatrix> multiplyByMatrix(std::shared_ptr<aMatrix> other) override;

    virtual std::shared_ptr<aMatrix> transpose() const override;

    virtual std::shared_ptr<aVector> multiplyByVector(std::shared_ptr<aVector> vector) override;

    virtual void softCopy(std::shared_ptr<aMatrix> other) override;

    virtual void hardCopy(std::shared_ptr<aMatrix> other) override;

    virtual bool isZero() const override;

    virtual void dump(std::string filename) const override;

    virtual BlockMatrix* clone() const override;

    void resize(const unsigned int& nRows, const unsigned int& nCols);

    std::shared_ptr<aMatrix> block(const unsigned int& iblock, const unsigned int& jblock) const override;

    std::shared_ptr<BlockMatrix> getSubmatrix(const unsigned int& ibegin, const unsigned int& iend,
                                 const unsigned int& jbegin, const unsigned int& jend) const;

    void setBlock(const unsigned int& iblock, const unsigned int& jblock,
                  std::shared_ptr<aMatrix> matrix) override;

    // inline bool isFinalized() const {return M_isFinalized;}

    void printPattern() const;

    void close() override;

    void open() override {M_isOpen = true;}

    inline bool isOpen() const {return M_isOpen;}

    inline void checkOpen() const {if (!isOpen()) throw new Exception("BlockMatrix must be open for this operation!");}

    inline void checkClosed() const {if (isOpen()) throw new Exception("BlockMatrix must be closed for this operation!");}

    // void finalize() {};
    //
    // std::shared_ptr<aMatrix> collapse() const {};

protected:
    BlockMatrix();

    void updateNormInf();
    // I introduce this only because on mac the operator+= behaves weirdly
    // void sumMatrix(const BlockMatrix& other) {};
    //
    // void multiplyCoeff(const double& coeff) {};

    Grid            M_matrixGrid;
    bool            M_isOpen;
};

}

#endif // BLOCKMATRIX_HPP
