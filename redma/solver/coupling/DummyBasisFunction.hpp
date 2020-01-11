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

// these are taken from https://arxiv.org/pdf/1701.02709.pdf (formula 3.8)

#ifndef DUMMYBASISFUNCTION_HPP
#define DUMMYBASISFUNCTION_HPP

#include <functional>

#include <redma/solver/coupling/BasisFunctionFunctor.hpp>
#include <redma/utils/Exception.hpp>

namespace RedMA
{

class DummyBasisFunction : public BasisFunctionFunctor
{
public:
    DummyBasisFunction(const GeometricFace& face,
                       std::string type);

    return_Type operator()(const Vector3D& pos) override;
};

}  // namespace RedMA

#endif  // DUMMYBASISFUNCTION_HPP
