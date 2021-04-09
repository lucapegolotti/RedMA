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

#ifndef SEGMENTATIONSMERGER_HPP
#define SEGMENTATIONSMERGER_HPP

#include <redma/RedMA.hpp>
#include <redma/geometry/SegmentationParser.hpp>
#include <redma/geometry/TreeStructure.hpp>

#include <lifev/core/filter/GetPot.hpp>

namespace RedMA
{

/*! \brief Merger of multiple SegmentationParsers.
 *
 * At the moment, we support only two segmentations.
 */
class SegmentationsMerger
{
    typedef LifeV::VectorSmall<3>                      Vector3D;
    typedef GeometricFace                              Contour;
    typedef LifeV::MatrixSmall<3,3>                    Matrix3D;
    typedef shp<SegmentationParser>                    SegmentationParserPtr;

public:

    /*! \brief Constructor.
     *
     * \param The DataContainer of the problem.
     * \param comm The MPI Communicator.
     * \param verbose If true, output is printed to standard output.
     */
    SegmentationsMerger(const DataContainer& datafile,
                        EPETRACOMM comm,
                        bool verbose);

    /*! \brief Merge multiple SegmentationParsers.
     *
     * \param parsers The shared pointers SegmentationParsers.
     * \return A TreeStructure.
     */
    TreeStructure merge(std::vector<SegmentationParserPtr> parsers);

    /*! \brief Merge two segmentations.
     *
     * \param segmentationFather Father segmentation.
     * \param segmentationChild Child segmentation.
     * \param outputTree Resulting tree.
     */
    void mergeTwoSegmentations(SegmentationParserPtr segmentationFather,
                               SegmentationParserPtr segmentationChild,
                               TreeStructure& outputTree);

    /*! \brief Returns connectivity matrix of the segmentations.
     *
     * The connectivity matrix is 1 in (i,j) if the ith segmentation is the
     * father of the jth one.
     *
     * \return The connectivity matrix.
     */
    bool** getConnectivity();

private:
    double findClosestPoint(const Vector3D& targetPoint,
                            SegmentationParserPtr toSearchIn,
                            unsigned int startIndexSearch,
                            unsigned int& indexOfClosestPoint);

    double placeBifurcation(Vector3D initialCenter,
                            Vector3D initialAxis,
                            Vector3D initialTransverse,
                            double initialRadius,
                            SegmentationParserPtr segmentationFather,
                            SegmentationParserPtr segmentationChild,
                            shp<BifurcationSymmetric> bifurcation);

    double rotateBifurcation(Vector3D initialCenter,
                             Vector3D initialAxis,
                             Vector3D initialTransverse,
                             double initialRadius,
                             SegmentationParserPtr segmentationFather,
                             SegmentationParserPtr segmentationChild,
                             shp<BifurcationSymmetric> bifurcation);

    // params: bx,by,bz,rotation_axis_x,rotation_axis_y,rotation_axis_z,
    // alpha, scale, alpha_axis, out1_alphax, out1_alphay, out1_alphaz,
    // out2_alphax, out2_alphay, out2_alphaz
    int deformBifurcation(shp<BifurcationSymmetric> bifurcation,
                          std::vector<double> params);

    // params: alpha_axis, out1_alphax, out1_alphay, out1_alphaz,
    // out2_alphax, out2_alphay, out2_alphaz
    int deformPlacedBifurcation(shp<BifurcationSymmetric> bifurcation,
                                std::vector<double> params);

    double computeLoss(shp<BifurcationSymmetric> bifurcation,
                       SegmentationParserPtr segmentationFather,
                       SegmentationParserPtr segmentationChild,
                       const double inletConst, const double outletConst);

    double optimizeLoss(shp<BifurcationSymmetric> bifurcation,
                        SegmentationParserPtr segmentationFather,
                        SegmentationParserPtr segmentationChild,
                        std::vector<double>& params,
                        double lambda,
                        const double tol, const unsigned int nMax);

    void deallocateConnectivity(bool** connectivity, unsigned int size);

    DataContainer       M_datafile;
    EPETRACOMM          M_comm;
    bool                M_verbose;
    double              M_constCenters;
    double              M_constNormals;
    double              M_constRadius;
};

}  // namespace RedMA

#endif  // SEGMENTATIONSMERGER_HPP
