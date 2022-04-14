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

#ifndef SNAPSHOTSSTEADYSAMPLER_HPP
#define SNAPSHOTSSTEADYSAMPLER_HPP

#include <redma/RedMA.hpp>
#include <redma/problem/DataContainer.hpp>
#include <redma/problem/GlobalProblem.hpp>
#include <redma/reduced_basis/LatinHypercube.hpp>
#include <redma/geometry/GeometryPrinter.hpp>

#include <cmath>
#include <iomanip>
#include <fstream>
#include <random>

namespace RedMA
{

/// Class handling the generation of the steady snapshots.
class SnapshotsSteadySampler
{
public:
    SnapshotsSteadySampler(const DataContainer& data, EPETRACOMM comm, unsigned int numSamples);

    /// Take the snapshots.
    void takeSnapshots(const unsigned int& Nstart = 0);

    inline void setInflow(const std::function<double(double,double,double)>& inflow) {M_inflow=inflow;};

    void dumpSnapshots(GlobalProblem& problem, std::string outdir, const std::vector<double>& params);

    void printCurrentSample(std::map<std::string, double> sample);

    std::vector<double> getParametersValuesAsVector(const std::map<std::string, double>& sample);

    std::map<std::string, double> getCurrentSample(unsigned int i);

private:
    LatinHypercube                                      M_LHS;
    DataContainer                                       M_data;
    EPETRACOMM                                          M_comm;
    std::function<double(double,double,double)>         M_inflow;
};

}  // namespace RedMA

#endif  // SNAPSHOTSSTEADYSAMPLER_HPP
