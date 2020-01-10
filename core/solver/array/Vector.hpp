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

#ifndef VECTOR_HPP
#define VECTOR_HPP

#include <Exception.hpp>
#include <AbstractVector.hpp>

namespace RedMA
{

template<class InVectorType>
class Vector : public AbstractVector
{
    typedef std::shared_ptr<InVectorType>            InVectorTypePtr;
public:
    Vector();

    // soft copy
    Vector(const InVectorTypePtr& vector);

    // soft copy
    Vector(const Vector<InVectorType>& other);

    // hard copy
    Vector<InVectorType>& operator=(const Vector<InVectorType>& other);

    // hard copy
    Vector<InVectorType>& operator=(const InVectorTypePtr& vector);

    Vector<InVectorType>& operator+=(const Vector<InVectorType>& other);

    Vector<InVectorType>& operator*=(const double& coeff);

    InVectorTypePtr& get();

    InVectorTypePtr get() const;

    bool isZero() const;

    bool isNull() const;

private:
    InVectorTypePtr             M_inVector;
};

}  // namespace RedMA

#include <Vector_imp.hpp>

#endif  // VECTOR_HPP
