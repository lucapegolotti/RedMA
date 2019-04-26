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

#ifndef ROSENBROCKALGORITHM_HPP
#define ROSENBROCKALGORITHM_HPP

#include <TimeMarchingAlgorithm.hpp>

namespace RedMA
{

template <class AssemblerType>
class RosenbrockAlgorithm : public TimeMarchingAlgorithm<AssemblerType>
{
protected:
    typedef GlobalAssembler<AssemblerType>      GlobalAssemblerType;

public:
    RosenbrockAlgorithm(const GetPot& datafile);

    virtual void solveTimestep(const double &time, double &dt,
                               const GlobalAssemblerType& assembler);

};

}  // namespace RedMA

#include <RosenbrockAlgorithm_imp.hpp>

#endif  // ROSENBROCKALGORITHM_HPP
