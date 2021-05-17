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

#include <redma/array/aVector.hpp>
#include <redma/array/aMatrix.hpp>

#include <fstream>

namespace RedMA
{

class aFunctionProvider
{
public:

    aFunctionProvider() {};

    virtual shp<aVector> getZeroVector() const = 0;

    virtual shp<aMatrix> getMass(const double& time, const shp<aVector>& sol) = 0;

    virtual shp<aMatrix> getPressureMass(const double& time, const shp<aVector>& sol) = 0;

    virtual shp<aMatrix> getMassJacobian(const double& time, const shp<aVector>& sol) = 0;

    virtual shp<aVector> getRightHandSide(const double& time, const shp<aVector>& sol) = 0;

    virtual shp<aMatrix> getJacobianRightHandSide(const double& time, const shp<aVector>& sol) = 0;

    virtual void apply0DirichletBCs(shp<aVector> vector) const = 0;

    virtual void applyDirichletBCs(const double& time, shp<aVector> vector) const = 0;

    virtual void setExtrapolatedSolution(const shp<aVector>& exSol) = 0;
};

}

#endif // aFUNCTIONPROVIDER_HPP
