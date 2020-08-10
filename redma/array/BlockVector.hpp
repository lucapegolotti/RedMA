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

#ifndef BLOCKVECTOR_HPP
#define BLOCKVECTOR_HPP

#include <redma/utils/Exception.hpp>
#include <redma/array/DenseVector.hpp>

#include <boost/numeric/ublas/matrix.hpp>

#include <math.h>

namespace RedMA
{

class BlockVector : public aVector
{
    typedef boost::numeric::ublas::matrix<std::shared_ptr<aVector>>        Grid;

public:

    BlockVector();

    BlockVector(const BlockVector& other);

    BlockVector(const unsigned int& nRows);

    virtual void add(std::shared_ptr<aVector> other) override;

    virtual void multiplyByScalar(const double& coeff) override;

    virtual void dump(std::string namefile) const override;

    virtual void softCopy(std::shared_ptr<aVector> other) override;

    virtual void hardCopy(std::shared_ptr<aVector> other) override;

    virtual aVector* clone() const override;

    virtual bool isZero() const override;

    std::string getString(const char& delimiter) const override;

    double norm2() const override;

    void setBlock(const unsigned int& iblock, std::shared_ptr<aVector> vector);

    // BlockVector operator*(const double& coeff) const;
    //
    // BlockVector& operator*=(const double& coeff);
    //
    // BlockVector& operator+=(const BlockVector& other);
    //
    // BlockVector& operator-=(const BlockVector& other);
    //
    // BlockVector operator+(const BlockVector& other) const;
    //
    // BlockVector operator-(const BlockVector& other) const;

    std::shared_ptr<aVector> block(const unsigned int& iblock) const;

    std::shared_ptr<BlockVector> getSubvector(const unsigned int& ibegin, const unsigned int& iend) const;

    void resize(const unsigned int& nRows);

    void updateNormInf();

    inline void close() {M_isOpen = false;}

    inline void open() {M_isOpen = true;}

    inline bool isOpen() const {return M_isOpen;}

    inline void checkOpen() const {if (!isOpen()) throw new Exception("BlockMatrix must be open for this operation!");}

    inline void checkClosed() const {if (isOpen()) throw new Exception("BlockMatrix must be closed for this operation!");}

protected:
    Grid            M_vectorGrid;
    bool            M_isOpen;
};

}

#endif // BLOCKVECTOR_HPP
