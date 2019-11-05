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

#include <sstream>

#define THRESHOLD 3e-2

namespace RedMA
{

// we can use the geometric face for the contour as well
class SegmentationParser
{
    typedef LifeV::VectorSmall<3>                          Vector3D;
    typedef GeometricFace                                  Contour;
    typedef std::shared_ptr<Epetra_Comm>                   commPtr_Type;
    typedef LifeV::MatrixSmall<3,3>                        Matrix3D;

public:
    // interpolation method specifies how to associate contour information
    // to path points in between two contours
    SegmentationParser(commPtr_Type comm, std::string pthName,
                       std::string ctgrName, std::string interpolationMethod,
                       bool verbose);

    ~SegmentationParser();

    void traversePath(std::string pthName);

    void traverseSegmentation(std::string ctgrName);

    TreeStructure createTree(const int& lengthTubes,
                             const double& constVector,
                             const double& constNormal,
                             int indexBegin = -1, int indexEnd = -1);

private:
    Vector3D get3DVectorFromXMLElement(tinyxml2::XMLElement* data);

    void linearInterpolation();

    // alpha is angle about axis, theta is bending angle
    // The function that maps the outlet center c is the following
    // c_tilde = rot(alpha) scale * A bend(c,theta) + b,
    // where rot is the rotation matrix about the axis, A2 and b come from the
    // affine transformation, and bend is the bend function. L is the length
    // of the tube to be deformed in the default configuration.
    // Const1 and const2 are used in the loss function
    double Fbending(const double& alpha, const double& theta,
                    const Matrix3D& A, const Vector3D& b,
                    const double& scale,
                    const double& L, const Vector3D& axis,
                    const Contour& target,
                    const double& const1,
                    const double& const2);

    double aPosterioriCheck(const double& alpha, const double& theta,
                            const Matrix3D& A, const Vector3D& b,
                            const double& scale,
                            const double& L, const Vector3D& axis,
                            const Contour& target);

    double optimizeBending(double& alpha, double& theta,
                           const double& maxIt, const double& tol,
                           const Matrix3D& A, const Vector3D& b,
                           const double& scale,
                           const double& L, const Vector3D& axis,
                           const Contour& target,
                           const double& const1,
                           const double& const2);

    void bend(const double& alpha, const double& theta,
              Vector3D& center, Vector3D& normal,
              const Matrix3D& A, const Vector3D& b,
              const double& scale,
              const double& L, const Vector3D& axis);

    std::vector<Contour>    M_contours;
    std::vector<Vector3D>   M_path;
    std::vector<Vector3D>   M_tangents; // corresponding to ith path point
    std::vector<Vector3D>   M_rotation; // corresponding to ith path point
    std::vector<Contour>    M_contoursComplete;
    bool                    M_verbose;
    commPtr_Type            M_comm;
    unsigned int            M_indexBegin;
    unsigned int            M_indexEnd;
    std::vector<double>     M_cumulativeDistance;
};

}  // namespace RedMA

#endif  // SEGMENTATIONPARSER_HPP
