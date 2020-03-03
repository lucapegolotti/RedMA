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
#include <redma/solver/problem/DataContainer.hpp>
#include <redma/solver/problem/ProblemFEM.hpp>

#include <redma/geometry/GeometryPrinter.hpp>

#include <boost/filesystem.hpp>

namespace RedMA
{

class SnapshotsSampler
{
public:
    SnapshotsSampler(const DataContainer& data, EPETRACOMM comm);

    void takeSnapshots();

    void dumpSnapshots(ProblemFEM& problem, std::string outdir);

private:
    DataContainer       M_data;
    EPETRACOMM          M_comm;
};

}  // namespace RedMA

#endif  // SNAPSHOTSSAMPLER_HPP