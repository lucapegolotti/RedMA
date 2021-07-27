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

#ifndef GEOMETRYPRINTER_HPP
#define GEOMETRYPRINTER_HPP

#include <stdio.h>
#include <string>
#include <queue>

#include <tinyxml2.h>

#include <redma/geometry/TreeStructure.hpp>

namespace RedMA
{

/// Class to print a TreeStructure to a .xml file.
class GeometryPrinter
{
public:

    /// Default constructor.
    GeometryPrinter();

    /*! \brief Save a TreeStructure to a file.
     *
     * \param The TreeStructure.
     * \param name Name of the file.
     * \param comm The MPI Communicator.
     */
    void saveToFile(TreeStructure& tree,
                    std::string name,
                    shp<Epetra_Comm> comm);
};

}  // namespace RedMA

#endif  // GEOMETRYPRINTER_HPP
