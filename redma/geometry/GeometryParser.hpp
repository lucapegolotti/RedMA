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

#include <redma/utils/PrintLog.hpp>
#include <redma/geometry/BuildingBlock.hpp>
#include <redma/geometry/Tube.hpp>
#include <redma/geometry/BifurcationSymmetric.hpp>
#include <redma/geometry/Aorta.hpp>
#include <redma/geometry/AortaBifurcation0.hpp>
#include <redma/geometry/AortaBifurcation1.hpp>
#include <redma/geometry/TreeStructure.hpp>

#include <Epetra_SerialComm.h>
#include <Epetra_MpiComm.h>

namespace RedMA
{

class GeometryParser
{
    typedef shp<Epetra_Comm>  commPtr_Type;
    typedef shp<BuildingBlock>  BuildingBlockPtr;
    typedef tinyxml2::XMLElement  XMLEl;

public:
    GeometryParser(const GetPot& datfile, std::string fileName,
                   commPtr_Type comm, bool verbose);

    void traverseXML(XMLEl* curElement,
                     unsigned int IDfather);

    TreeStructure& getTree();

private:
    BuildingBlockPtr parseElement(const XMLEl* element, unsigned int& outletParent);

    commPtr_Type M_comm;
    TreeStructure M_tree;
    GetPot  M_datafile;
    bool M_verbose;
    int M_maxNumBlocks;
    unsigned int M_numBlocks;
};

}  // namespace RedMA

#endif  // GEOMETRYPARSER_HPP
