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
#include <redma/array/DistributedVector.hpp>
#include <redma/array/Wrap.hpp>

#include <boost/numeric/ublas/matrix.hpp>

#include <math.h>

namespace RedMA
{

class BlockVector : public aVector
{
    typedef boost::numeric::ublas::matrix<shp<aVector>>        Grid;

public:

    BlockVector();

    virtual ~BlockVector() {};

    BlockVector(const unsigned int& nRows);

    virtual void add(shp<aVector> other) override;

    virtual void multiplyByScalar(const double& coeff) override;

    virtual void dump(std::string namefile) const override;

    virtual void shallowCopy(shp<aDataWrapper> other) override;

    virtual void deepCopy(shp<aDataWrapper> other) override;

    virtual BlockVector* clone() const override;

    virtual bool isZero() const override;

    std::string getString(const char& delimiter) const override;

    double norm2() const override;

    void setBlock(const unsigned int& iblock, shp<aVector> vector);

    bool globalTypeIs(Datatype type);

    unsigned int level();

    shp<BlockVector> convertInnerTo(Datatype type, shp<Epetra_Comm> comm = nullptr);

    shp<aVector> block(const unsigned int& iblock) const override;

    shp<BlockVector> getSubvector(const unsigned int& ibegin, const unsigned int& iend) const;

    void resize(const unsigned int& nRows);

    void findComm();

    shp<Epetra_Comm> commPtr() {findComm(); return M_comm;}

    void copyPattern(shp<BlockVector> other, bool verbose = true);

    virtual double operator()(unsigned int index) override {throw new Exception("operator() undefined for BlockVector");}

    virtual shp<void> data() const override {return nullptr;};

    virtual void setData(shp<void>) override {throw new Exception("setData undefined for BlockVector");};

    virtual Datatype type() const override {return BLOCK;}

protected:
    shp<Epetra_Comm>              M_comm;
    Grid                          M_vectorGrid;
};

}

#endif // BLOCKVECTOR_HPP
