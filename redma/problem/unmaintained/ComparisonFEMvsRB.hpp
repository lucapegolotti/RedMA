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

#ifndef COMPARISONFEMVSRB_HPP
#define COMPARISONFEMVSRB_HPP

#include <redma/problem/GlobalProblem.hpp>
#include <redma/problem/ProblemRB.hpp>

namespace RedMA
{

class ComparisonFEMvsRB
{
public:
    ComparisonFEMvsRB(const DataContainer& data, EPETRACOMM comm);

    void setGlobalProblem(const shp<GlobalProblem>& problemFEM) {M_problemFEM = problemFEM;}

    void setProblemRB(const shp<ProblemRB>& problemRB) {M_problemRB = problemRB;}

    void runFEM();

    void runRB();

    void exportError();

    void dumpFEMSolution(std::string outdir);

    void loadFEMSolution(std::string indir);

    void exportFEM(unsigned int saveEvery);

    void exportRB(unsigned int saveEvery);

    double getTimeFem() {return M_timeFEM;}

    double getTimeRB() {return M_timeRB;}

private:
    shp<GlobalProblem>                                     M_problemFEM;
    shp<ProblemRB>                                      M_problemRB;
    std::vector<BlockVector>                            M_loadedSolutions;
    double                                              M_timeFEM;
    double                                              M_timeRB;
    DataContainer                                       M_data;
    EPETRACOMM                                          M_comm;
};

}

#endif // COMPARISONFEMVSRB_HPP
