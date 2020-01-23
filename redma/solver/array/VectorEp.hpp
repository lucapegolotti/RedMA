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

#ifndef VECTOREP_HPP
#define VECTOREP_HPP

#include <redma/RedMA.hpp>
#include <redma/solver/array/aMatrix.hpp>

namespace RedMA
{

class VectorEp : public aMatrix
{
public:
    VectorEp();

    VectorEp operator+(const VectorEp& other);

    VectorEp& operator+=(const VectorEp& other);

    VectorEp& operator*=(const double& coeff);

    void hardCopy(const VectorEp& other);

private:
    std::shared_ptr<VECTOREPETRA>  M_vector;
};

}

#endif // VECTOREP_HPP