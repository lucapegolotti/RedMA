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

#ifndef RBBASESMANAGER_HPP
#define RBBASESMANAGER_HPP

#include <redma/reduced_basis/RBBases.hpp>

#include <redma/RedMA.hpp>

namespace RedMA
{

/// Manager class for multiple reduced bases.
class RBBasesManager
{
public:
    /*! \brief Constructor.
     *
     * \param dataContainer A DataContainer.
     * \param comm A MPI Communicator.
     * \param idmeshmap A std map with key = ID (of building block), value = meshname.
     */
    RBBasesManager(const DataContainer& dataContainer,
                   EPETRACOMM comm,
                   std::map<unsigned int, std::string> idmeshmap);

    /// Load the singular values from file
    void loadSingularValues();

    /// Load the reduced bases from file.
    void loadBases();

    /*! \brief Get reduced bases associated with a mesh.
     *
     * \param dataContainer A DataContainer.
     * \param comm A MPI Communicator.
     * \param idmeshmap A std map with key = ID (of building block), value = meshname.
     */
    shp<RBBases> getRBBases(std::string meshName);

private:
    EPETRACOMM                                      M_comm;
    DataContainer                                   M_data;
    std::map<unsigned int, std::string>             M_IDMeshMap;
    std::map<std::string, shp<RBBases>>             M_bases;
};

}  // namespace RedMA

#endif  // RBBASESMANAGER_HPP
