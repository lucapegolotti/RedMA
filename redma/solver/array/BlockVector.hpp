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
#include <redma/solver/array/BlockMatrix.hpp>

#include <boost/numeric/ublas/matrix.hpp>

namespace RedMA
{

template <class InVectorType>
class BlockVector : public BlockMatrix<InVectorType>
{
public:
    BlockVector();

    BlockVector(const unsigned int& nRows);

    InVectorType& block(const unsigned int& iblock);

    InVectorType block(const unsigned int& iblock) const;
};

}

#include "BlockVector_imp.hpp"

#endif // BLOCKVECTOR_HPP
