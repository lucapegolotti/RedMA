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

#ifndef BUILDINGBLOCK_HPP
#define BUILDINGBLOCK_HPP

#include <map>
#include <memory>
#include <sys/stat.h>
#include <sys/types.h>
#include <iostream>

#include <redma/utils/Exception.hpp>
#include <redma/utils/PrintLog.hpp>
#include <redma/geometry/GeometricParametersHandler.hpp>
#include <redma/geometry/MembraneThicknessComputer.hpp>
#include <redma/problem/DataContainer.hpp>

#include <lifev/core/mesh/RegionMesh.hpp>
#include <lifev/core/mesh/MeshData.hpp>
#include <lifev/core/mesh/MeshPartitioner.hpp>
#include <lifev/core/filter/GetPot.hpp>

#include <lifev/core/LifeV.hpp>
#include <lifev/core/array/VectorEpetra.hpp>
#include <Epetra_SerialComm.h>
#include <Epetra_MpiComm.h>
#include <lifev/core/mesh/MeshUtility.hpp>
#include <lifev/core/array/VectorSmall.hpp>
#include <lifev/core/array/MatrixSmall.hpp>
#include <lifev/core/filter/ExporterHDF5.hpp>
#include <lifev/core/filter/ExporterVTK.hpp>
#include <lifev/core/fem/FESpace.hpp>

#include <redma/RedMA.hpp>

namespace RedMA
{

/// A geometric face.
class GeometricFace
{
public:
    typedef LifeV::VectorSmall<3>   Vector3D;

    /// Default constructor.
    GeometricFace();

    /*! \brief Constructor.
     *
     * \param center Center of the face.
     * \param normal Normal of the face.
     * \param radius The radius of the face.
     * \param flag The flag of the face.
     */
    GeometricFace(Vector3D center,
                  Vector3D normal,
                  double radius,
                  unsigned int flag);

    /*! \brief Constructor.
     *
     * \param center Center of the face.
     * \param normal Normal of the face.
     * \param radius The radius of the face.
     * \param flag The flag of the face.
     * \param flagDisk The flag of the ring.
     */
    GeometricFace(Vector3D center,
                  Vector3D normal,
                  double radius,
                  unsigned int flag,
                  unsigned int flagRing);

    /// Print information regarding the geometric face.
    void print() const;

    /*! \brief Equality operator between geometric faces.
     *
     * \param lhs Left Hand Side face
     * \param rhs Right Hand Side face
     */
    friend bool operator==(const GeometricFace& lhs, const GeometricFace& rhs);

    unsigned int    M_flag;
    unsigned int    M_ringFlag;
    Vector3D        M_center;
    Vector3D        M_normal;
    double          M_radius;
};

/// Abstract building block containing the mesh of a subdomain.
class BuildingBlock
{
protected:

    typedef LifeV::VectorSmall<3>                           Vector3D;
    typedef LifeV::MatrixSmall<3,3>                         Matrix3D;
    typedef LifeV::ExporterVTK<MESH>                        Exporter;
    typedef LifeV::MeshUtility::MeshTransformer<MESH>       Transformer;
    typedef shp<GeometricParameter>                         GeometricParameterPtr;

public:

    /*! \brief Constructor.
     *
     * \param comm The MPI Communicator.
     * \param refinement The refinement level.
     * \param verbose If true, output is pushed to standard output.
     */
    BuildingBlock(EPETRACOMM comm,
                  std::string refinement,
                  bool verbose);

    /*! \brief Set the value of a parameter.
     *
     * \param key The key of the parameter.
     * \param value The value of the parameter.
     */
    int setParameterValue(std::string key,
                          double value);

    /*! \brief Get a parameter.
     *
     * \param name The name of the parameter.
     * \return Shared pointer to the geometric parameter.
     */
    GeometricParameterPtr getParameter(std::string name) {return M_parametersHandler.getParameter(name);};

    /*! \brief Get the map of all parameters.
     *
     * \return Parameters map (key = name, value = pointer to the GeometricParameter)
     */
    std::map<std::string,shp<GeometricParameter>>& getParametersMap();

    /*! \brief Getter for the GeometricParametersHandler.
     *
     * \return Reference to the GeometricParameterHandler.
     */
    GeometricParametersHandler& getGeometricParametersHandler() {return M_parametersHandler;}

    /*! \brief Read the mesh.
     *
     * \param meshdir The mesh where the meshes are stored.
     */
    void readMesh(std::string meshdir = "../../../meshes/");

    /*! \brief Number of children the building block accepts.
     *
     * \return The expected number of children.
     */
    virtual unsigned int expectedNumberOfChildren() = 0;

    /*! \brief Get the name of the BuildingBlock.
     *
     * \return The name.
     */
    std::string name();

    /*! \brief Apply the affine transformation.
     *
     * Apply a scaling, rotation, and translation to the mesh.
     *
     * \param transformMesh If true, the mesh is modified and transformed.
     */
    void applyAffineTransformation(bool transformMesh = true);

    /*! \brief Apply nonaffine transformation.
     *
     * This function is abstract.
     *
     * \param transformMesh If true, the mesh is modified and transformed.
     */
    virtual void applyNonAffineTransformation(bool transformMesh = true) = 0;

    /*! \brief Apply non-affine and affine transformation.
     *
     * \param transformMesh If true, the mesh is modified and transformed.
     */
    void applyGlobalTransformation(bool transformMesh = true);

    /*! \brief Dump the mesh to file.
     *
     * \param outdir The output directory.
     * \param meshdir The directory where the meshes are stored.
     * \param outputName Name of the mesh to be dumped.
     */
    void dumpMesh(std::string outdir,
                  std::string meshdir,
                  std::string outputName);

    /*! \brief Get an outlet.
     *
     * \param indexFace The index of the face.
     * \return The outlet geometric face.
     */
    GeometricFace getOutlet(unsigned int indexFace) const;

    /*! \brief Get the outlets.
     *
     * \return Vector of the outlet geometric faces.
     */
    std::vector<GeometricFace> getOutlets() const;

    /*! \brief Get an inlet.
     *
     * \param indexInlet The index of the inlet.
     * \return The inlet geometric face.
     */
    GeometricFace getInlet(unsigned int indexFace) const;

    /*! \brief Get the inlets.
     *
     * \return Vector of the inlet geometric faces.
     */
    std::vector<GeometricFace> getInlets() const;

    /*! \brief Set the affine transformation such that the inlet matches a face.
     *
     * \param parentOutlet The face to match.
     */
    void mapChildInletToParentOutlet(GeometricFace parentOutlet);

    /*! \brief Set boolean to specify if this block is a child.
     *
     * \param The boolean to set.
     */
    void setIsChild(bool isChild);

    /*! \brief Getter for the mesh.
     *
     * \retrun Shared pointer to the mesh.
     */
    shp<MESH> getMesh() const;

    /*! \brief Set the datafile.
     *
     * \param The datafile to set.
     */
    void setDatafile(const DataContainer& datafile);

    /// Randomize parameters.
    void setRandom();

    /// Set the inlet and outlets.
    virtual void resetInletOutlets() = 0;

    /*! \brief Compute the rotation and angle to reach a vector.
     *
     * \param vectorToMove The starting vector.
     * \param vectorToReach The final vector.
     * \param axis The rotation axis (output).
     * \param alphax The rotation angle (output).
     */
    static void computeRotationAxisAndAngle(Vector3D vectorToMove,
                                            Vector3D vectorToReach,
                                            Vector3D& axis,
                                            double& alphax);

    /*! \brief Compute a rotation matrix.
     *
     * \param axis The rotation axis index (0 = x, 1 = y, 2 = z)
     * \param angle The rotation angle.
     * \return The rotation matrix.
     */
    static Matrix3D computeRotationMatrix(unsigned int axis,
                                          double angle);

    /*! \brief Compute a rotation matrix.
     *
     * \param axis The rotation axis.
     * \param angle The rotation angle.
     * \return The rotation matrix.
     */
    static Matrix3D computeRotationMatrix(Vector3D axis,
                                          double angle);

    /// Compute the thickness of the membrane over the building block lateral surface.
    void computeMembraneThickness();

    /*! \brief Getter of the membrane thickness.
     *
     * \return A vector storing the thickness of the membrane over the building block lateral surface.
     */
    inline shp<VECTOREPETRA> getMembraneThickness() const {return M_membraneThicknessComputer->getThickness();}

    /*! \brief Get the refinement.
     *
     * \param The refinement.
     */
    std::string getRefinement(){return M_refinement;};

    /*! \brief Get building blocks dependent parameters.
     *
     * \param index Index of the parameter.
     * \return The parameter name.
     */
    virtual std::string getOptionalParameter(unsigned int index){return "";};

    /*! \brief Get MPI Communicator.
     *
     * \return The MPI Communicator.
     */
    inline EPETRACOMM getComm() const {return M_comm;}

    /*! \brief Returns true if the building block is child.
     *
     * \return True if the building block is child.
     */
    inline bool getIsChild() const {return M_isChild;}

    /*! \brief Get the mesh name.
     *
     * \return The mesh name.
     */
    inline std::string getMeshName() const {return M_meshName;}

    /*! \brief Get the discretization method.
     *
     * \return The discretization method.
     */
    inline std::string getDiscretizationMethod() const {return M_discrMethod;}

    /*! \brief Get the assembler type.
     *
     * \return The assembler type.
     */
    inline std::string getAssemblerType() const {return M_assemblerType;}

    /*! \brief Getter of the wall flag.
     *
     * \return The wall flag.
     */
    inline unsigned int getWallFlag() {return M_wallFlag;}

    /*! \brief Set the discretization method.
     *
     * \param method The discretization method.
     */
    inline void setDiscretizationMethod(std::string method) {M_discrMethod = method;}

    /*! \brief Set the assembler type.
     *
     * \param method The assembler type.
     */
    inline void setAssemblerType(std::string assembler) {M_assemblerType = assembler;}

    /*! \brief Compute the Jacobian global transformation.
     *
     * \param x First component of the point in which the Jacobian must be computed.
     * \param y Second component of the point in which the Jacobian must be computed.
     * \param z Third component of the point in which the Jacobian must be computed.
     * \return The Jacobian matrix.
     */
    Matrix3D computeJacobianGlobalTransformation(const double& x,
                                                 const double& y,
                                                 const double& z);

    /*! \brief Compute the Jacobian non affine transformation.
     *
     * \param x First component of the point in which the Jacobian must be computed.
     * \param y Second component of the point in which the Jacobian must be computed.
     * \param z Third component of the point in which the Jacobian must be computed.
     * \return The Jacobian matrix.
     */
    virtual Matrix3D computeJacobianNonAffineTransformation(const double& x,
                                                            const double& y,
                                                            const double& z) = 0;

    /*! \brief Compute the inverse of a 3x3 matrix.
     *
     * \param matrix Input matrix
     * \return The inverse.
     */
    static void matrixInverse(Matrix3D& matrix,
                              Matrix3D* inverse);

    /*! \brief Apply global transformation to 3D vector.
     *
     * \param x In/out parameter of x-component.
     * \param y In/out parameter of y-component.
     * \param z In/out parameter of z-component.
     */
    void globalTransf(double& x,
                      double& y,
                      double& z);

    /*! \brief Apply affine transformation to 3D vector.
     *
     * \param x In/out parameter of x-component.
     * \param y In/out parameter of y-component.
     * \param z In/out parameter of z-component.
     */
    void affineTransf(double& x,
                      double& y,
                      double& z);

    /*! \brief Apply nonaffine transformation to 3D vector.
     *
     * \param x In/out parameter of x-component.
     * \param y In/out parameter of y-component.
     * \param z In/out parameter of z-component.
     */
    virtual void nonAffineTransf(double& x,
                                 double& y,
                                 double& z) {};

    /*! \brief Compute a rotation matrix with respect to a local (n, t1, t2) reference system
     *
     * \param vec 3D vector, representing the normal vector
     * \param index Index of the normal component in the resulting local reference system It defaults to 2.
     * \return The rotation matrix to the local reference system
     */
    Matrix3D computeLocalRotationMatrix(Vector3D vec, const unsigned int& index = 2) const;

    /*! \brief Compute two orthonormal tangent vectors to a given input vector
     *
     * \param vec 3D vector of reference
     * \return A pair of tangent orthonormal vectors to the input vector
     */
    static std::pair<Vector3D, Vector3D> computeLocalTangentVersors(const Vector3D& normal);

    shp<VECTOREPETRA> getDisplacement() { return M_displacement; };

    void setDefaultStiffness(shp<MATRIXEPETRA> stiff);

    shp<MATRIXEPETRA> getDefaultStiffness() { return M_originalStiffness; };

protected:
    void applyAffineTransformationGeometricFace(GeometricFace& face,
                                                const Matrix3D& affineMatrix,
                                                const Vector3D& translation,
                                                const double& scale);

    static void rotationFunction(double& x, double& y, double& z,
                                 const Matrix3D& affMatrix,
                                 const Vector3D& transl, const double& scale);

    static double fZero(const double& t, const double& x, const double& y,
                        const double& z, const LifeV::ID& i);

    GeometricParametersHandler            M_parametersHandler;
    std::string                           M_name;

    std::string                           M_refinement;
    std::string                           M_meshName;
    shp<MESH>                             M_mesh;

    Matrix3D                              M_R;
    Matrix3D                              M_Raxis;

    EPETRACOMM                            M_comm;

    bool                                  M_verbose;

    std::string                           M_datafileName;

    std::vector<GeometricFace>            M_inlets;
    std::vector<GeometricFace>            M_outlets;

    bool                                  M_isChild;

    double                                M_inletScale;
    double                                M_scale;

    Vector3D                              M_inletTranslation;
    Vector3D                              M_inletRotationAxis;
    Vector3D                              M_translation;
    double                                M_inletAngle;

    DataContainer                         M_datafile;

    std::string                           M_discrMethod;
    std::string                           M_assemblerType;

    unsigned int                          M_wallFlag;

    shp<MATRIXEPETRA>                     M_originalStiffness;

    shp<MembraneThicknessComputer>        M_membraneThicknessComputer;

    shp<VECTOREPETRA>                     M_displacement;
};

}  // namespace RedMA

#endif  // BUILDINGBLOCK_HPP
