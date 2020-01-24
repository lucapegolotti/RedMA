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

#ifndef LINEARSYSTEMSOLVER_HPP
#define LINEARSYSTEMSOLVER_HPP

#include <redma/RedMA.hpp>
#include <redma/solver/array/BlockVector.hpp>
#include <redma/solver/array/BlockMatrix.hpp>

#include <memory>

namespace RedMA
{

template<class InVectorType, class InMatrixType>
class LinearSystemSolver
{
    typedef BlockVector<InVectorType>               BV;
    typedef BlockMatrix<InMatrixType>               BM;

public:
    LinearSystemSolver(const GetPot& datafile);

    BV solve(BM matrix, BV rh);

private:
    GetPot                  M_datafile;
};

}

#include "LinearSystemSolver_imp.hpp"

#endif // LINEARSYSTEMSOLVER_HPP
