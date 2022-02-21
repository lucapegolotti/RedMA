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

#ifndef aPROBLEM_HPP
#define aPROBLEM_HPP

#include <redma/problem/DataContainer.hpp>

namespace RedMA
{

/// Generic interface for a mathematical problem.
class aProblem
{
public:
    /*! \brief Constructor.
     *
     * \param data Datafile.
     */
    aProblem(const DataContainer& data);

    /// Virtual setup method
    virtual void setup() = 0;

    /// Virtual solve method
    virtual void solve() = 0;

    /// Virtual solve method for the steady problem
    virtual void solveSteady() = 0;

    /// Getter for the datafile
    DataContainer& getData() {return M_data;}

protected:
    DataContainer  M_data;
};

}

#endif // aPROBLEM_HPP
