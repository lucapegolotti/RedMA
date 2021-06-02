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

#include <redma/problem/aProblem.hpp>
#include <redma/array/BlockVector.hpp>
#include <redma/geometry/TreeStructure.hpp>
#include <redma/geometry/GeometryParser.hpp>

#include <redma/solver/time_marching_algorithms/TimeMarchingAlgorithmFactory.hpp>
#include <redma/assemblers/block/BlockAssembler.hpp>
#include <redma/assemblers/DefaultAssemblersLibrary.hpp>

#include <memory>

namespace RedMA
{

/*! \brief Object managing a problem defined on multiple subdomains.
 *
 * This class provides an interface to solve a differential problem on multiple
 * subdomains.
 */
class GlobalProblem : public aProblem
{
    typedef shp<BlockVector>                   BV;
    typedef shp<BlockMatrix>                   BM;
public:
    /*! \brief Constructor.
     *
     * Here, the geometry is parsed from the .xml contained in the datafile
     * and, if doSetup == true, the setup method is launched.
     *
     * \param data Datafile with problem settings.
     * \param comm MPI Communicator
     * \param doSetup If true, setup is launched.
     */
    GlobalProblem(const DataContainer& data, EPETRACOMM comm, bool doSetup = true);

    /*! \brief Setup method.
     *
     * Read reference meshes and deform them, set default assemblers (these are used,
     * e.g, in the Piola transformation), initialize time marching algorithm.
     */
    virtual void setup();

    /*! \brief Solve method.
     *
     * Method that launches the time loop over the timesteps.
     */
    virtual void solve();

    /*! \brief Exporting method
     *
     * Method to export solutions already stored in .txt files
     */
    void exportFromFiles(const std::string& path);

    /*! \brief Set if solutions should be stored at each timestep.
     *
     * The solutions can be retrieved with getSolutions().
     */
    void doStoreSolutions() {M_storeSolutions = true;}

    /*! \brief Getter for the geometric tree.
     *
     * \return Reference to M_tree.
     */
    inline TreeStructure& getTree() {return M_tree;}

    /*! \brief Get stored solutions.
     *
     *  \return Vector of shared pointer to BlockVectors, containing the solutions.
     */
    inline std::vector<BV> getSolutions() {return M_solutions;}

    /*! \brief Get a vector with all the timesteps.
     *
     * \return A vector with all the timesteps.
     */
    inline std::vector<double> getTimesteps() {return M_timestepsSolutions;}

    /*! \brief Getter for the BlockAssembler.
     *
     * \return Shared pointer to the internal BlockAssembler.
     */
    inline shp<BlockAssembler> getBlockAssembler() {return M_assembler;}

    /*! \brief Check if internal BlockAssembler is composed of only FE assemblers.
     *
     * \return True if all the primal assemblers are finite element assemblers.
     */
    bool isFEProblem();

    /*! \brief Check if internal BlockAssembler is composed of only RB assemblers.
     *
     * \return True if all the primal assemblers are reduced basis assemblers.
     */
    bool isRBProblem();

private:

    shp<aTimeMarchingAlgorithm>                M_TMAlgorithm;
    shp<BlockAssembler>                        M_assembler;
    BV                                         M_solution;
    TreeStructure                              M_tree;
    GeometryParser                             M_geometryParser;
    bool                                       M_storeSolutions;
    std::vector<BV>                            M_solutions;
    std::vector<double>                        M_timestepsSolutions;
    shp<DefaultAssemblersLibrary>              M_defaultAssemblers;
    EPETRACOMM                                 M_comm;
};

}

#endif // PROBLEMFEM_HPP
