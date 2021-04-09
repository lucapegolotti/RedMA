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

#ifndef GEOMETRYPARSER_HPP
#define GEOMETRYPARSER_HPP

#include <stdio.h>
#include <string>
#include <tinyxml2.h>

#include <redma/RedMA.hpp>

#include <redma/geometry/building_blocks/BuildingBlock.hpp>
#include <redma/geometry/building_blocks/Tube.hpp>
#include <redma/geometry/building_blocks/BifurcationSymmetric.hpp>
#include <redma/geometry/building_blocks/Aorta.hpp>
#include <redma/geometry/building_blocks/AortaBifurcation0.hpp>
#include <redma/geometry/building_blocks/AortaBifurcation1.hpp>
#include <redma/geometry/TreeStructure.hpp>

namespace RedMA
{

/// Generator of a TreeStructure given an input file.
class GeometryParser
{
    typedef shp<BuildingBlock>  BuildingBlockPtr;
    typedef tinyxml2::XMLElement  XMLEl;

public:

    /*! \brief Constructor.
     *
     * \param datafile A datafile.
     * \param fileName The geometry file.
     * \param comm The MPI Communicator
     */
    GeometryParser(const DataContainer& datafile,
                   std::string fileName,
                   EPETRACOMM comm,
                   bool verbose);

    /*! \brief Traverse XML file.
     *
     * \param curElement Current XMLElement
     *
     */
    void traverseXML(XMLEl* curElement,
                     unsigned int IDfather);

    /*! \brief Return reference to the TreeStructure.
     *
     * \return Reference to the TreeStructure.
     */
    TreeStructure& getTree();

private:
    BuildingBlockPtr parseElement(const XMLEl* element,
                                  unsigned int& outletParent);

    EPETRACOMM          M_comm;
    TreeStructure       M_tree;
    DataContainer       M_datafile;
    bool                M_verbose;
    int                 M_maxNumBlocks;
    unsigned int        M_numBlocks;
};

}  // namespace RedMA

#endif  // GEOMETRYPARSER_HPP
