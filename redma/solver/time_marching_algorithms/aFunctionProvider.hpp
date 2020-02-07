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

#ifndef aFUNCTIONPROVIDER_HPP
#define aFUNCTIONPROVIDER_HPP

#include <redma/RedMA.hpp>

#include <redma/solver/array/BlockVector.hpp>
#include <redma/solver/array/BlockMatrix.hpp>

#include <fstream>

namespace RedMA
{

template <class InVectorType, class InMatrixType>
class aFunctionProvider
{
public:

    aFunctionProvider() {};

    virtual BlockVector<InVectorType> getZeroVector() const = 0;

    virtual BlockMatrix<InMatrixType> getMass(const double& time,
                                      const BlockVector<InVectorType>& sol) = 0;

    virtual BlockMatrix<InMatrixType> getMassJacobian(const double& time,
                                      const BlockVector<InVectorType>& sol) = 0;

    virtual BlockVector<InVectorType> getRightHandSide(const double& time,
                                      const BlockVector<InVectorType>& sol) = 0;

    virtual BlockMatrix<InMatrixType> getJacobianRightHandSide(const double& time,
                                      const BlockVector<InVectorType>& sol) = 0;

};

}

#endif // aFUNCTIONPROVIDER_HPP
