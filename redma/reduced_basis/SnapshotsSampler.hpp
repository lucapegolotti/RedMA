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

#ifndef SNAPSHOTSSAMPLER_HPP
#define SNAPSHOTSSAMPLER_HPP

#include <redma/RedMA.hpp>
#include <redma/problem/DataContainer.hpp>
#include <redma/problem/GlobalProblem.hpp>

#include <redma/geometry/GeometryPrinter.hpp>

#include <cmath>
#include <iomanip>
#include <fstream>
#include <random>

namespace RedMA
{

/// Class handling the generation of the snapshots.
class SnapshotsSampler
{
public:
    SnapshotsSampler(const DataContainer& data, EPETRACOMM comm);

    /// Take the snapshots.
    void takeSnapshots(const unsigned int& Nstart = 0);

    inline void setInflow(const std::function<double(double,double,double)>& inflow) {M_inflow=inflow;};

    void dumpSnapshots(GlobalProblem& problem, std::string outdir, const std::vector<double> array_params);

    void transformSnapshotsWithPiola(std::string snapshotsDir,
                                     unsigned int fieldIndex,
                                     unsigned int maxSnapshot);

    std::vector<double> inflowSnapshots(const std::vector<std::vector<double>>& param_bounds);

private:
    DataContainer                                       M_data;
    EPETRACOMM                                          M_comm;

    std::function<double(double,double,double)>         M_inflow;
};

}  // namespace RedMA

#endif  // SNAPSHOTSSAMPLER_HPP
