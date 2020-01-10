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

#include <Exception.hpp>

#include <Vector.hpp>
#include <lifev/core/array/VectorEpetra.hpp>

namespace RedMA
{

class BlockVector : public AbstractVector
{
    typedef LifeV::VectorEpetra                                  VectorEpetra;
    typedef std::shared_ptr<VectorEpetra>                        VectorEpetraPtr;
public:
    BlockVector();

    BlockVector(unsigned int numRows);

    void resize(unsigned int numRows);

    Vector<VectorEpetra>& block(unsigned int index);

    Vector<VectorEpetra> block(unsigned int index) const;

    VectorEpetraPtr& operator[](unsigned int index);

    void push_back(Vector<VectorEpetra> newElement);

    unsigned int size() const;

private:
    std::vector<Vector<VectorEpetra> >    M_vectors;
};

}  // namespace RedMA

#endif  // BLOCKVECTOR_HPP
