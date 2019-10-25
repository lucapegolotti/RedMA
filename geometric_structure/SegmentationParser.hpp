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

#ifndef SEGMENTATIONPARSER_HPP
#define SEGMENTATIONPARSER_HPP

#include <BuildingBlock.hpp>
#include <TreeStructure.hpp>
#include <Tube.hpp>

#include <Epetra_SerialComm.h>
#include <Epetra_MpiComm.h>

#include <tinyxml2.h>

namespace RedMA
{

// we can use the geometric face for the contour as well
class SegmentationParser
{
    typedef LifeV::VectorSmall<3>                          Vector3D;
    typedef GeometricFace                                  Contour;
    typedef std::shared_ptr<Epetra_Comm>                   commPtr_Type;


public:
    SegmentationParser(commPtr_Type comm, std::string pthName,
                       std::string ctgrName, bool verbose);

    ~SegmentationParser();

    void traversePath(std::string pthName);

    void traverseSegmentation(std::string ctgrName);

    TreeStructure createTree(int indexBegin = -1, int indexEnd = -1);

private:
    Vector3D get3DVectorFromXMLElement(tinyxml2::XMLElement* data);

    std::vector<Contour>    M_contours;
    std::vector<Vector3D>   M_path;
    std::vector<Vector3D>   M_tangents; // corresponding to ith path point
    std::vector<Vector3D>   M_rotation; // corresponding to ith path point
    bool                    M_verbose;
    commPtr_Type            M_comm;
};

}  // namespace RedMA

#endif  // SEGMENTATIONPARSER_HPP
