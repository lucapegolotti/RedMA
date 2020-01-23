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

#ifndef PROBLEMFEM_HPP
#define PROBLEMFEM_HPP

#include <redma/RedMA.hpp>

#include <redma/solver/problem/aProblem.hpp>
#include <redma/solver/array/BlockVector.hpp>

#include <redma/solver/time_marching_algorithms/TimeMarchingAlgorithmFactory.hpp>
#include <redma/solver/assemblers/AssemblerFactory.hpp>

#include <memory>

namespace RedMA
{

class ProblemFEM : public aProblem
{
public:
    ProblemFEM(const GetPot& datafile);

    virtual void setup();

    virtual void solve();

    void solveTimestep(const double& t, double& dt);

private:
    SHP(aTimeMarchingAlgorithm<BlockVector<FEVECTOR> >)  M_timeMarchingAlgorithm;
    SHP(aAssembler)              M_assembler;
    BlockVector<FEVECTOR>        M_solution;
};

}

#endif // PROBLEMFEM_HPP
