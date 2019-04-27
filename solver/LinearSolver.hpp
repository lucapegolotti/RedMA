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

#ifndef LINEARSOLVER_HPP
#define LINEARSOLVER_HPP

#include <lifev/core/array/VectorEpetraStructured.hpp>
#include <lifev/core/array/MatrixEpetraStructured.hpp>

namespace RedMA
{

class LinearSolver
{
protected:
    typedef LifeV::VectorEpetraStructured                Vector;
    typedef std::shared_ptr<Vector>                      VectorPtr;
    typedef LifeV::MatrixEpetraStructured<double>        Matrix;
    typedef std::shared_ptr<Matrix>                      MatrixPtr;

public:
    LinearSolver();

    unsigned int solve(VectorPtr& solution, MatrixPtr matrix, VectorPtr rhs) const;

};

}  // namespace RedMA

#endif  // LINEARSOLVER_HPP
