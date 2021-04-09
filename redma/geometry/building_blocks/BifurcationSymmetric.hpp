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

#ifndef BIFURCATIONSYMMETRIC_HPP
#define BIFURCATIONSYMMETRIC_HPP

#include <redma/geometry/building_blocks/BuildingBlock.hpp>
#include <redma/geometry/building_blocks/NonAffineDeformer.hpp>

namespace RedMA
{

/// Building block of a symmetric bifurcation.
class BifurcationSymmetric : public BuildingBlock
{
public:

    /*! \brief Default constructor.
     *
     * \param comm The MPI Communicator.
     * \param refinement The refinement.
     * \param verbose If true, output is pushed to standard output.
     * \param angle The angle of the bifurcation.
     * \param randomizible If true, the geometrical parameters are randomizible.
     */
    BifurcationSymmetric(EPETRACOMM comm,
                         std::string refinement = "coarse",
                         bool verbose = false,
                         int angle = 50,
                         bool randomizable = true);

     /*! \brief Return the expected number of children.
      *
      * \return The expected number of children (2).
      */
    virtual inline unsigned int expectedNumberOfChildren() override
    {
        return 2;
    }

    /*! \brief Apply nonaffine transformation.
     *
     * This functions performs a rotation of the outlets by solving a linear
     * elasticity problem in which we impose the outlets displacement as boundary
     * conditions.
     *
     * \param transformMesh If true, the mesh is modified and transformed.
     */
    virtual void applyNonAffineTransformation(bool transformMesh = true) override;

    /*! \brief Get the center for the rotation of the outlets.
     *
     * The center is based on the outlets position.
     *
     * \return The center.
     */
    inline Vector3D getCenter() const {return M_center;};

    /*! \brief Get normal to the inlet.
     *
     * \return The normal to the inlet.
     */
    inline Vector3D getInletNormal() const {return M_inletNormalRef;};

    /*! \brief Get normal to the plane in which the bifurcation layes.
     *
     * \return The desired vector.
     */
    inline Vector3D getTransverse() const {return M_transverse;};

    /*! \brief Get the radius of the inlet.
     *
     * \return The radius of the inlet.
     */
    inline double getInletRadius() const {return M_inletRadiusRef;};

    /// Set the inlet and outlets.
    void resetInletOutlets() override;

    /*! \brief Get building blocks dependent parameters.
     *
     * If index = 0, returns the angle of the bifurcation.
     *
     * \param index Index of the parameter.
     * \return The parameter name.
     */
    virtual std::string getOptionalParameter(unsigned int index) override;

    /*! \brief Compute the Jacobian non affine transformation.
     *
     * We approximate this Jacobian with the identity.
     *
     * \param x First component of the point in which the Jacobian must be computed.
     * \param y Second component of the point in which the Jacobian must be computed.
     * \param z Third component of the point in which the Jacobian must be computed.
     * \return The Jacobian matrix.
     */
    virtual Matrix3D computeJacobianNonAffineTransformation(const double& x,
                                                            const double& y,
                                                            const double& z) override {return M_identity3D;};

private:
    static double outletMapFunction(const double& t, const double& x,
                               const double& y, const double& z,
                               const LifeV::ID& i,
                               const GeometricFace& targetFace,
                               const Vector3D& desiredCenter,
                               const Matrix3D& rotationMatrix);

    void bend(const double& out1_alphax,
              const double& out1_alphay,
              const double& out1_alphaz,
              const double& out2_alphax,
              const double& out2_alphay,
              const double& out2_alphaz,
              shp<Transformer> transformer,
              bool transformMesh = true);

    void rotateGeometricFace(const GeometricFace& face, Vector3D& rotatedCenter,
                             Vector3D& rotatedNormal,
                             const Matrix3D& rotationMatrix,
                             const Vector3D& rotationCenter);

    void computeCenter();

    std::string getStringMesh(std::string refinement);

    Vector3D M_inletCenterRef;
    Vector3D M_inletNormalRef;
    Vector3D M_outlet1CenterRef;
    Vector3D M_outlet1NormalRef;
    Vector3D M_outlet2CenterRef;
    Vector3D M_outlet2NormalRef;

    Vector3D M_center;

    Vector3D M_transverse;

    Matrix3D M_identity3D;

    double M_inletRadiusRef;
    double M_outlet1RadiusRef;
    double M_outlet2RadiusRef;
    int    M_angle;
};

}  // namespace RedMA

#endif  // BIFURCATIONSYMMETRIC_HPP
