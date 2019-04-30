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

#ifndef TIMEMARCHINGALGORITHMFACTORY_HPP
#define TIMEMARCHINGALGORITHMFACTORY_HPP

#include <memory>

#include <TimeMarchingAlgorithm.hpp>
#include <RosenbrockAlgorithm.hpp>
#include <Exception.hpp>

namespace RedMA
{

template <class AssemblerType>
std::shared_ptr<TimeMarchingAlgorithm<AssemblerType> >
TimeMarchingAlgorithmsFactory(const GetPot& datafile,
                              GlobalAssembler<AssemblerType>* assembler)
{
    std::string marchingAlgorithmString =
        datafile("time_discretization/algorithm", "rosenbrock");

    if (!std::strcmp(marchingAlgorithmString.c_str(), "rosenbrock"))
    {
        typedef RosenbrockAlgorithm<AssemblerType>  ReturnType;
        std::shared_ptr<ReturnType>
                returnPtr(new RosenbrockAlgorithm<AssemblerType>(datafile,
                                                                 assembler));

        return returnPtr;
    }
    else
    {
        std::string errorMsg = "Time marching algorithm of type " +
                    marchingAlgorithmString + " is not implemented!";

        throw Exception(errorMsg);
    }
    return nullptr;
}


}  // namespace RedMA

#endif  // TIMEMARCHINGALGORITHMFACTORY_HPP
