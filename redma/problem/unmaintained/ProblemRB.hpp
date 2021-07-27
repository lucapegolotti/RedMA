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

#ifndef PROBLEMRB_HPP
#define PROBLEMRB_HPP

#include <redma/RedMA.hpp>

#include <redma/problem/aProblem.hpp>
#include <redma/array/BlockVector.hpp>
#include <redma/geometry/TreeStructure.hpp>
#include <redma/geometry/GeometryParser.hpp>

#include <redma/solver/time_marching_algorithms/TimeMarchingAlgorithmFactory.hpp>
#include <redma/assemblers/block/BlockAssembler.hpp>
#include <redma/solver/assemblers/DefaultAssemblersLibrary.hpp>

#include <memory>

namespace RedMA
{

class ProblemRB : public aProblem
{
    typedef shp<BlockVector>                   BV;
    typedef shp<BlockVector>                   BBV;
    typedef shp<BlockMatrix>                   BM;
    typedef shp<BlockMatrix>                   BBM;
public:
    ProblemRB(const DataContainer& data, EPETRACOMM comm, bool doSetup = true);

    virtual void setup();

    virtual void solve();

    void doStoreSolutions(){M_storeSolutions = true;}

    void doStoreNonLinearTerms(){M_storeNonLinearTerms = true;}

    inline TreeStructure& getTree() {return M_tree;}

    inline std::vector<BBV> getSolutions() {return M_solutions;}

    inline std::vector<BBV> getNonLinearTerms() {return M_nonLinearTerms;}

    inline std::vector<double> getTimesteps() {return M_timestepsSolutions;}

    inline shp<BlockAssembler> getBlockAssembler() {return M_assembler;}

private:
    shp<aTimeMarchingAlgorithm>               M_TMAlgorithm;
    shp<BlockAssembler>                       M_assembler;
    BBV                                       M_solution;
    TreeStructure                             M_tree;
    GeometryParser                            M_geometryParser;
    bool                                      M_storeSolutions;
    bool                                      M_storeNonLinearTerms;
    std::vector<BBV>                          M_solutions;
    std::vector<BBV>                          M_nonLinearTerms;
    std::vector<double>                       M_timestepsSolutions;
};

}

#endif // PROBLEMRB_HPP
