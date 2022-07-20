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

#ifndef BYPASS_HPP
#define BYPASS_HPP

#include <redma/geometry/building_blocks/BuildingBlock.hpp>
#include <redma/geometry/building_blocks/NonAffineDeformer.hpp>

namespace RedMA
{

/// Building block of a femoropopliteal bypass.
    class Bypass : public BuildingBlock
    {
    public:
        /*! \brief Default constructor.
         *
         * \param comm The MPI Communicator.
         * \param name The name of the mesh.
         * \param verbose If true, output is pushed to standard output.
         * \param boundary_layer If true, the mesh with the boundary layer is considered
         * \param randomizable If true, the geometrical parameters are randomizable.
         */
        Bypass(EPETRACOMM comm,
               std::string name = "bypass",
               bool verbose = false,
               bool boundary_layer = true,
               bool isBifurcation = true,
               unsigned int activeStenosis = 1,
               bool randomizable = true
               );

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

        /// Set the active stenosis.
        /// i = 0 activates the old one
        /// i = 1 activates the one on the direction of outlet zero
        /// i = 2 activates the one on the direction of outlet one
        /// i = 3 activates the one close to the inlet
        void setActiveStenosis(unsigned int i);

        /*! \brief Compute the Jacobian non affine transformation.
         *
         * We approximate it with the identity.
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

        static double inletMapFunction(const double& t, const double& x,
                                       const double& y, const double& z,
                                       const LifeV::ID& i,
                                       const GeometricFace& targetFace,
                                       const Vector3D& desiredCenter,
                                       const Matrix3D& rotationMatrix);

        void bend(const double& in1_alphax,
                  const double& in1_alphay,
                  const double& in1_alphaz,
                  const double& in2_alphax,
                  const double& in2_alphay,
                  const double& in2_alphaz,
                  shp<Transformer> transformer,
                  bool transformMesh = true);

        void rotateGeometricFace(const GeometricFace& face, Vector3D& rotatedCenter,
                                 Vector3D& rotatedNormal,
                                 const Matrix3D& rotationMatrix,
                                 const Vector3D& rotationCenter);

        static double stenosisDeformation(const double &z);

        static double stenosisBC(const double &t, const double &x,
                                 const double &y, const double &z,
                                 const LifeV::ID &i, const double& amplitude, const double& width,
                                 const Vector3D& stenosisCentre, const Vector3D& stenosisNormal,
                                 const Matrix3D& distorsionMatrix);

        void addStenosis(const double &amplitude, const double &width,
                    shp <Transformer> transformer, bool transformMesh);

        void computeCenter();

        void computeStenosisCenter();

        void computeStenosisOuterNormal();

        void setStenosisAttributes();

        std::map<unsigned int, std::map<std::string, Vector3D>> M_stenosisAttributes;

        Vector3D M_inletCenterRef1;
        Vector3D M_inletNormalRef1;
        Vector3D M_inletCenterRef2;
        Vector3D M_inletNormalRef2;
        Vector3D M_outletCenterRef;
        Vector3D M_outletNormalRef;

        Vector3D M_center;
        Vector3D M_stenosisCenter;
        Vector3D M_stenosisOuterNormal;

        Vector3D M_Eigenvector1;
        Vector3D M_Eigenvector2;
        Vector3D M_Eigenvector3;

        double M_diameterAtStenosis;

        bool M_isBifurcation;

        void setDistorsionMatrix();

        Matrix3D M_identity3D;
        Matrix3D M_distorsionMatrix;

        double M_inletRadiusRef1;
        double M_inletRadiusRef2;
        double M_outletRadiusRef;
    };

}  // namespace RedMA

#endif  // BYPASS_HPP
