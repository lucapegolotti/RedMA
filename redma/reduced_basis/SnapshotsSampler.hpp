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
#include <functional>
#include <redma/RedMA.hpp>
#include <redma/problem/DataContainer.hpp>
#include <redma/problem/GlobalProblem.hpp>

#include <redma/geometry/GeometryPrinter.hpp>

namespace RedMA
{
    typedef std::function<double(double)> myFunctionType;
/// Class handling the generation of the snapshots.
class SnapshotsSampler
{
public:
    /*! \brief Constructor.
     *
     * \param dataContainer A DataContainer.
     * \param comm A MPI Communicator.
     */
    SnapshotsSampler(const DataContainer& data,
                     EPETRACOMM comm,bool geo);
    SnapshotsSampler(const DataContainer& data,
                     EPETRACOMM comm,bool geo,myFunctionType inflow_1,myFunctionType inflow_2, double index_, double index_max_);


    /// Take the snapshots.
    void takeSnapshots();

    /*! \brief Dump snapshots to file.
     *
     * \param problem A GlobalProblem containing the solutions.
     * \param outdir Output directory
     */
    void dumpSnapshots(GlobalProblem& problem,
                       std::string outdir);

private:
    DataContainer       M_data;
    EPETRACOMM          M_comm;
    bool                geometry;
    bool                analyticInflow;
    std::function<double(double)> inflow1;
    std::function<double(double)> inflow2;
    double     index;
    double     index_max;
};

}  // namespace RedMA

#endif  // SNAPSHOTSSAMPLER_HPP
