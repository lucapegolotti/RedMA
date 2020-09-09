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

#ifndef aVECTOR_HPP
#define aVECTOR_HPP

#include <redma/utils/Exception.hpp>

#include <memory>

#define ZEROTHRESHOLD   1e-15

namespace RedMA
{

enum Datatype{DENSE, SPARSE, DISTRIBUTED, BLOCK};

class aVector
{
public:
    aVector(Datatype type) :
      M_nRows(0),
      M_type(type),
      M_normInf(0.0)
    {
    }

    virtual ~aVector() {};

    virtual void add(std::shared_ptr<aVector> other) = 0;

    virtual void multiplyByScalar(const double& coeff) = 0;

    virtual void dump(std::string namefile) const = 0;

    virtual void softCopy(std::shared_ptr<aVector> other) = 0;

    virtual void hardCopy(std::shared_ptr<aVector> other) = 0;

    virtual aVector* clone() const = 0;

    virtual bool isZero() const = 0;

    virtual std::string getString(const char& delimiter) const = 0;

    virtual double norm2() const = 0;

    inline double normInf() const {return M_normInf;}

    inline unsigned int nRows() const {return M_nRows;}

    inline Datatype type() {return M_type;}

protected:
    template <class Type>
    void checkType(std::shared_ptr<Type> other, Datatype type)
    {
        if (other->type() != type)
            throw new Exception("Unexpected datatype!");
    }

    unsigned int        M_nRows;
    Datatype            M_type;
    double              M_normInf;
};

}

#endif // aVECTOR_HPP
