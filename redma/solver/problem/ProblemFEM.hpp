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
#include <redma/geometry/TreeStructure.hpp>
#include <redma/geometry/GeometryParser.hpp>

#include <redma/solver/time_marching_algorithms/TimeMarchingAlgorithmFactory.hpp>
#include <redma/solver/assemblers/BlockAssembler.hpp>

#include <memory>

namespace RedMA
{

class ProblemFEM : public aProblem
{
    typedef BlockVector<FEVECTOR>                   BV;
    typedef BlockVector<BV>                         BBV;
    typedef BlockMatrix<FEMATRIX>                   BM;
    typedef BlockMatrix<BM>                         BBM;
public:
    ProblemFEM(const DataContainer& data, EPETRACOMM comm);

    virtual void setup();

    virtual void solve();

    void solveTimestep(const double& t, double& dt);

private:
    SHP(aTimeMarchingAlgorithm<BV COMMA BM>)  M_TMAlgorithm;
    SHP(aAssembler<BV COMMA BM>)              M_assembler;
    BBV                                       M_solution;
    TreeStructure                             M_tree;
    GeometryParser                            M_geometryParser;
};

}

#endif // PROBLEMFEM_HPP
