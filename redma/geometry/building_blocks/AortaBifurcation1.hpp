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

#ifndef AORTABIFURCATION1_HPP
#define AORTABIFURCATION1_HPP

#include <redma/geometry/building_blocks/BuildingBlock.hpp>
#include <redma/geometry/building_blocks/NonAffineDeformer.hpp>

namespace RedMA
{

/// Building block of a bifurcation in the aorta.
class AortaBifurcation1 : public BuildingBlock
{
public:

    /*! \brief Default constructor.
     *
     * \param comm The MPI Communicator.
     * \param refinement The refinement.
     * \param name The name of the mesh.
     * \param verbose If true, output is pushed to standard output.
     */
    AortaBifurcation1(EPETRACOMM comm,
                      std::string refinement,
                      std::string name = "aortabif1",
                      bool verbose = false);

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
     * This function does not do anything.
     *
     * \param transformMesh If true, the mesh is modified and transformed.
     */
    virtual void applyNonAffineTransformation(bool transformMesh) override;

    /*! \brief Get building blocks dependent parameters.
     *
     * This function does not do anything.
     *
     * \param index Index of the parameter.
     * \return The parameter name.
     */
    std::string getOptionalParameter(unsigned int index) override;

    /// Set the inlet and outlets.
    void resetInletOutlets() override;

    /*! \brief Compute the Jacobian non affine transformation.
     *
     * This function does not do anything.
     *
     * \param x First component of the point in which the Jacobian must be computed.
     * \param y Second component of the point in which the Jacobian must be computed.
     * \param z Third component of the point in which the Jacobian must be computed.
     * \return The Jacobian matrix.
     */
    virtual Matrix3D computeJacobianNonAffineTransformation(const double& x,
                                                            const double& y,
                                                            const double& z) override {};

private:

    Vector3D M_inletCenterRef;
    Vector3D M_inletNormalRef;
    Vector3D M_outletCenterRef1;
    Vector3D M_outletNormalRef1;
    Vector3D M_outletCenterRef2;
    Vector3D M_outletNormalRef2;

    double M_inletRadiusRef;
    double M_outletRadiusRef1;
    double M_outletRadiusRef2;
};

}  // namespace RedMA

#endif  // AORTA_HPP
