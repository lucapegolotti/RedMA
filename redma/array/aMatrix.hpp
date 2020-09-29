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

#ifndef aMATRIX_HPP
#define aMATRIX_HPP

#include <redma/array/aVector.hpp>
#include <redma/utils/Exception.hpp>

#include <memory>

namespace RedMA
{

class aMatrix
{
public:

    aMatrix(Datatype type) :
      M_nRows(0),
      M_nCols(0),
      M_type(type),
      M_normInf(0.0)
    {
    }

    virtual ~aMatrix() {};

    // for every operation, the priority is to do it in place (on *this) if possible.
    // The copy of *this is left to the user if necessary
    virtual void add(std::shared_ptr<aMatrix> other) = 0;

    virtual void multiplyByScalar(const double& coeff) = 0;

    virtual std::shared_ptr<aVector> multiplyByVector(std::shared_ptr<aVector> vector) = 0;

    virtual std::shared_ptr<aMatrix> multiplyByMatrix(std::shared_ptr<aMatrix> other) = 0;

    virtual std::shared_ptr<aMatrix> transpose() const = 0;

    virtual void dump(std::string namefile) const = 0;

    virtual void softCopy(std::shared_ptr<aMatrix> other) = 0;

    virtual void hardCopy(std::shared_ptr<aMatrix> other) = 0;

    virtual aMatrix* clone() const = 0;

    virtual bool isZero() = 0;

    virtual void setBlock(const unsigned int& i,
                          const unsigned int& j,
                          std::shared_ptr<aMatrix> matrix) {throw new Exception("Method setBlock(uint,uint,SHP(aMatrix)) not overloaded");}

    virtual std::shared_ptr<aMatrix> block(const unsigned int& i,
                                           const unsigned int& j) const {throw new Exception("Method block(uint,uint) not overloaded");}

    virtual std::shared_ptr<void> data() const {throw new Exception("Method data() not overloaded");}

    virtual void setData(std::shared_ptr<void> data) {throw new Exception("Method setData() not overloaded");}

    inline double normInf() const {return M_normInf;}

    inline unsigned int nRows() const {return M_nRows;}

    inline unsigned int nCols() const {return M_nCols;}

    inline Datatype type() {return M_type;}

    virtual void close() {}

    virtual void open() {}

protected:
    template <class Type>
    void checkType(std::shared_ptr<Type> other, Datatype type)
    {
        if (other->type() != type)
            throw new Exception("Unexpected datatype of matrix!");
    }

    unsigned int        M_nRows;
    unsigned int        M_nCols;
    Datatype            M_type;
    double              M_normInf;
};

}

#include "aMatrix_imp.hpp"

#endif // aMATRIX_HPP
