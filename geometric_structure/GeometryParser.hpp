// Reduced Modeling of Arteries (ReMA)
// Copyright (C) 2019  Luca Pegolotti
//
// ReMA is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// ReMA is distributed in the hope that it will be useful,
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

#include <PrintLog.hpp>
#include <BuildingBlock.hpp>
#include <Tube.hpp>
#include <BifurcationSymmetric.hpp>
#include <TreeStructure.hpp>

#include <Epetra_SerialComm.h>
#include <Epetra_MpiComm.h>

namespace ReMA
{

class GeometryParser
{
    typedef std::shared_ptr<Epetra_Comm>            commPtr_Type;

public:
    GeometryParser(std::string fileName, commPtr_Type comm, bool verbose);

    void traverseXML(tinyxml2::XMLElement* curElement, unsigned int IDfather);

    TreeStructure& getTree();

private:
    std::shared_ptr<BuildingBlock> parseElement(const tinyxml2::XMLElement* element);

    commPtr_Type M_comm;

    bool M_verbose;

    TreeStructure M_tree;
};

}  // namespace ReMA

#endif  // GEOMETRYPARSER_HPP
