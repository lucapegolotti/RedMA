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

#ifndef TIMEMARCHINGALGORITHM_HPP
#define TIMEMARCHINGALGORITHM_HPP

#include <GlobalAssembler.hpp>

namespace RedMA
{

template <class AssemblerType>
class TimeMarchingAlgorithm
{
protected:
    typedef GlobalAssembler<AssemblerType>      GlobalAssemblerType;

public:
    TimeMarchingAlgorithm(const GetPot& datafile);

    virtual void solveTimestep(const double &time, double &dt,
                               const GlobalAssemblerType& assembler) = 0;

protected:
    GetPot  M_datafile;
};

}  // namespace RedMA

#include <TimeMarchingAlgorithm_imp.hpp>

#endif  // TIMEMARCHINGALGORITHM_HPP
