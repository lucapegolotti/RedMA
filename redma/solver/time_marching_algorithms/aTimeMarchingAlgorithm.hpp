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

#ifndef aTIMEMARCHINGALGORITHM_HPP
#define aTIMEMARCHINGALGORITHM_HPP

#include <redma/RedMA.hpp>
#include <redma/solver/assemblers/aAssembler.hpp>
#include <redma/solver/array/BlockVector.hpp>
#include <redma/solver/system_solver/SystemSolver.hpp>

#include <lifev/core/filter/GetPot.hpp>

namespace RedMA
{

template <class InVectorType, class InMatrixType>
class aTimeMarchingAlgorithm
{
public:
    aTimeMarchingAlgorithm(const GetPot& datafile);

    virtual BlockVector<InVectorType> advance(const double& time, double& dt,
                       SHP(aAssembler<InVectorType COMMA InMatrixType>) assembler,
                       int& status) = 0;

protected:
    GetPot                                      M_datafile;
    SystemSolver<InVectorType, InMatrixType>    M_systemSolver;
};

}

#include "aTimeMarchingAlgorithm_imp.hpp"

#endif // aTIMEMARCHINGALGORITHM_HPP
