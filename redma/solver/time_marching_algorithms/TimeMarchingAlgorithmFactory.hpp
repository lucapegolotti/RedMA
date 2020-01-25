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

#include <redma/RedMA.hpp>
#include <redma/solver/time_marching_algorithms/aTimeMarchingAlgorithm.hpp>
#include <redma/solver/time_marching_algorithms/BDF.hpp>
#include <redma/utils/Exception.hpp>

#include <memory>

namespace RedMA
{

template <class InVectorType, class InMatrixType>
SHP(aTimeMarchingAlgorithm<InVectorType COMMA InMatrixType>)
TimeMarchingAlgorithmFactory(const GetPot& datafile);

}

#include "TimeMarchingAlgorithmFactory_imp.hpp"

#endif // TIMEMARCHINGALGORITHMFACTORY_HPP
