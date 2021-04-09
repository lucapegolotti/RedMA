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

#ifndef TUBE_HPP
#define TUBE_HPP

#include <redma/geometry/building_blocks/BuildingBlock.hpp>
#include <redma/geometry/building_blocks/NonAffineDeformer.hpp>

namespace RedMA
{

/// Building block of a tube.
class Tube : public BuildingBlock
{
public:

    /*! \brief Default constructor.
     *
     * \param comm The MPI Communicator.
     * \param refinement The refinement.
     * \param verbose If true, output is pushed to standard output.
     * \param diameter The diameter of the tube.
     * \param length The length of the tube.
     * \param randomizible If true, the geometrical parameters are randomizible.
     */
    Tube(EPETRACOMM comm,
         std::string refinement = "coarse",
         bool verbose = false,
         int diameter = 1,
         int length = 1,
         bool randomizable = true);

    /*! \brief Return the expected number of children.
     *
     * \return The expected number of children (1).
     */
    virtual inline unsigned int expectedNumberOfChildren() override
    {
        return 1;
    }

    /*! \brief Apply nonaffine transformation.
     *
     * This functions performs a scaling of the outlet, length, and a bending
     * in one direction of the tube.
     *
     * \param transformMesh If true, the mesh is modified and transformed.
     */
    virtual void applyNonAffineTransformation(bool transformMesh) override;

    /*! \brief Get refinement.
     *
     * \return The refinement.
     */
    std::string getStringMesh(std::string refinement);

    /*! \brief Get the radius of the inlet.
     *
     * \return The radius of the inlet.
     */
    inline double getInletRadius() const {return M_inletRadiusRef;};

    /*! \brief Get the normal to the inlet.
     *
     * \return The normal to the inlet.
     */
    Vector3D getInletNormal() {return M_inletNormalRef;};

    /*! \brief Get the length in the reference configuration.
     *
     * \return The length.
     */
    double getDefLength() {return (M_outletCenterRef-M_inletCenterRef).norm();};

    /*! \brief Get building blocks dependent parameters.
     *
     * If index = 0, returns the default diameter of the tube.
     * If index = 1, returns the default lenght of the tube.
     *
     * \param index Index of the parameter.
     * \return The parameter name.
     */
    std::string getOptionalParameter(unsigned int index) override;

    /// Set the inlet and outlets.
    void resetInletOutlets() override;

    /*! \brief Get the diameter of the tube.
     *
     * \return The diameter of the tube.
     */
    inline unsigned int getDiameter() {return M_diameter;};

    /*! \brief Get the length of the tube.
     *
     * \return The length of the tube.
     */
    inline unsigned int getLength() {return M_length;};

    /*! \brief Apply bend function to an input point.
     *
     * \param x First component of the input vector.
     * \param y Second component of the input vector.
     * \param z Third component of the input vector.
     * \param bendAngle The bending angle.
     * \param L The current length of the tube.
     */
    static void bendFunctionAnalytic(double& x,
                                     double& y,
                                     double& z,
                                     const double& bendAngle,
                                     const double& L);

    /*! \brief Compute the Jacobian non affine transformation.
     *
     *
     * \param x First component of the point in which the Jacobian must be computed.
     * \param y Second component of the point in which the Jacobian must be computed.
     * \param z Third component of the point in which the Jacobian must be computed.
     * \return The Jacobian matrix.
     */
    virtual Matrix3D computeJacobianNonAffineTransformation(const double& x,
                                                            const double& y,
                                                            const double& z) override;

    /*! \brief Apply nonaffine transformation to 3D vector.
     *
     * \param x In/out parameter of x-component.
     * \param y In/out parameter of y-component.
     * \param z In/out parameter of z-component.
     */
    virtual void nonAffineTransf(double& x,
                                 double& y,
                                 double& z) override;

private:
    void nonAffineScaling(const double& lengthRatio,
                          const double& outRadiusRatio,
                          shp<Transformer> transformer);

    static void scalingFunction(double& x, double& y, double& z,
                                const double& lenghtRatio,
                                const double& outRadiusRatio,
                                const double& L);

    Matrix3D computeJacobianNonAffineScaling(const double& x, const double& y, const double& z,
                                             const double& lengthRatio,
                                             const double& outRadiusRatio,
                                             const double& L);

    static double bendFunction(const double& t, const double& x, const double& y,
                               const double& z, const LifeV::ID& i,
                               const Vector3D& rotationCenter,
                               const Matrix3D& rotationMatrix);


    Matrix3D computeJacobianBend(const double& x, const double& y, const double& z,
                                 const double& bendAngle, const double& L);

    void bend(const double& bendAngle, shp<Transformer> transformer,
              bool transformMesh = true);

    Vector3D M_inletCenterRef;
    Vector3D M_inletNormalRef;
    Vector3D M_outletCenterRef;
    Vector3D M_outletNormalRef;

    double M_inletRadiusRef;
    double M_outletRadiusRef;

    unsigned int M_diameter;
    unsigned int M_length;

    Matrix3D M_jacobianScaling;
    Matrix3D M_jacobianBending;
};

}  // namespace RedMA

#endif  // TUBE_HPP
