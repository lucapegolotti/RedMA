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

#ifndef INTERFACEASSEMBLER_HPP
#define INTERFACEASSEMBLER_HPP

#include <redma/RedMA.hpp>

#include <redma/assemblers/abstract/aAssembler.hpp>
#include <redma/array/BlockMatrix.hpp>
#include <redma/array/BlockVector.hpp>
#include <redma/array/SparseMatrix.hpp>
#include <redma/array/DistributedVector.hpp>
#include <redma/utils/Exception.hpp>
#include <redma/coupling_basis_functions/BasisFunctionFactory.hpp>

#include <redma/reduced_basis/RBBases.hpp>

#include <lifev/eta/expression/Integrate.hpp>
#include <lifev/core/fem/DOFInterface3Dto3D.hpp>

#define THRESHOLDSTAB       1e-15

namespace RedMA
{

/*! \brief Class representing an interface between two subdomains.
 *
 */
class Interface
{
    typedef aAssembler                                   AssemblerType;

public:
    /// Default (empty) constructor.
    Interface ();

    /*! \brief Constructor taking various arguments.
     *
     * The arguments of the constructor are simply stored as members of the class.
     * We call father the upstream node and child the downstream one.
     *
     * \param assemblerFather Shared pointer to aAssembler of father node.
     * \param indexFather Index in the tree of the father node.
     * \param assemblerChild Shared pointer to aAssembler of child node.
     * \param indexChild Index in the tree of the child node.
     * \param interfaceID Index of the interface.
     */
    Interface(shp<AssemblerType> assemblerFather, const int& indexFather,
              shp<AssemblerType> assemblerChild, const int& indexChild,
              const unsigned int& interfaceID);

    shp<AssemblerType>      M_assemblerFather;
    shp<AssemblerType>      M_assemblerChild;
    int                     M_indexFather;
    int                     M_indexChild;
    unsigned int            M_ID;
    unsigned int            M_indexOutlet;
    unsigned int            M_indexInlet;
    unsigned int            M_interfaceFlag;
};

/*! \brief Class for the assembly of coupling matrices.
 *
 */
class InterfaceAssembler
{
    typedef aAssembler                                   AssemblerType;
    /*typedef LifeV::DOFInterface3Dto3D                    InterfaceType;
    typedef std::shared_ptr<InterfaceType>               InterfacePtrType;*/

    typedef LifeV::VectorSmall<3>                        Vector3D;
    typedef LifeV::MatrixSmall<3,3>                      Matrix3D;

public:

    /*! \brief Constructor taking a DataContainer as argument.
     *
     * \param data The DataContainer of the problem.
     * \param addNoSlipBC True in no-slip BCs at the vessel wall are desired (default)
     */
    InterfaceAssembler(const DataContainer& data,
                       const bool& addNoSlipBC = true);

    /*! \brief Constructor taking a DataContainer and an Interface as arguments.
     *
     * \param data The DataContainer of the problem.
     * \param interface The interface.
     * \param addNoSlipBC True in no-slip BCs at the vessel wall are desired (default)
     */
    InterfaceAssembler(const DataContainer& data,
                       const Interface& interface,
                       const bool& addNoSlipBC = true);

    /// Default empty destructor.
    virtual ~InterfaceAssembler() {}

    /*! \brief Setup method.
     *
     * This method allocates the coupling matrices and constructs them through
     * buildCouplingMatrices.
     */
    void setup();

    /*! \brief Build all the coupling matrices.
     *
     * Moreover, the matrices are projected onto the RB space if needed.
     */
    void buildCouplingMatrices();

    /*! \brief Build coupling matrix on a specific face and from a specific
     *         assembler.
     *
     * \param assembler Shared pointer to aAssembler.
     * \param face The GeometricFace.
     * \param matrixT Shared pointer to BlockMatrix of the transpose of the
     *                coupling matrix.
     * \param matrix Shared pointer to BlockMatrix of the coupling matrix.
     * \param isFather True if the considered block is the father one (default)
     *
     */
    void buildCouplingMatrices(shp<AssemblerType> assembler,
                               const GeometricFace& face,
                               shp<BlockMatrix> matrixT,
                               shp<BlockMatrix> matrix,
                               const bool isFather = true);

    /// Build stabilization matrix. Currently not implemented.
    void buildStabilizationMatrix(shp<AssemblerType> assembler,
                                  const GeometricFace& face,
                                  shp<BlockMatrix> matrix);

    /*! \brief Build Epetra map for the Lagrange multiplier space [Not implemented at the moment!].
     *
     * \param bfs Shared pointer to the BasisFunctionFunctor.
     */
    void buildMapLagrange(shp<BasisFunctionFunctor> bfs);

    /*! \brief Add coupling contribution to a right-hand side.
     *
     * \param time Current time (needed in derived classes).
     * \param rhs Shared pointer to the right-hand side.
     * \param sol Shared pointer to the solution.
     * \param nPrimalBlocks number of primal blocks in the problem.
     */
    virtual void addContributionRhs(const double& time,
                                    shp<BlockVector> rhs,
                                    shp<BlockVector> sol,
                                    const unsigned int& nPrimalBlocks);

    /*! \brief Add coupling contribution to the jacobian right-hand side.
     *
     * \param time Current time (needed in derived classes).
     * \param jac Shared pointer to the jacobian of the right-hand side.
     * \param sol Shared pointer to the solution.
     * \param nPrimalBlocks Number of primal blocks in the problem.
     */
    virtual void addContributionJacobianRhs(const double& time,
                                            shp<BlockMatrix> jac,
                                            shp<BlockVector> sol,
                                            const unsigned int& nPrimalBlocks);

    /*! \brief Getter for the interface.
     *
     * \return The Interface.
     */
    inline Interface getInterface() const {return M_interface;};

    /*! \brief Getter for the zero-vector.
     *
     * \return Shared pointer to BlockVector containing the zero vector.
     */
    shp<BlockVector> getZeroVector() const;

    /*! \brief Check norm of the stabilization term.
     *
     * \param Shared pointer to the solution.
     * \param nPrimalBlocks Number of primal blocks in the problem.
     */
    double checkStabilizationTerm(const shp<BlockVector>& sol,
                                  const unsigned int& nPrimalBlocks);

    /*! \brief Getter of the father's coupling matrix transpose.
     *
     * \return Father's coupling matrix transpose.
     */
    inline shp<BlockMatrix> getFatherBT() const {return M_fatherBT;}

    /*! \brief Getter of the father's coupling matrix.
     *
     * \return Father's coupling matrix.
     */
    inline shp<BlockMatrix> getFatherB() const {return M_fatherB;}

    /*! \brief Getter of the child's coupling matrix transpose.
     *
     * \return Child's coupling matrix transpose.
     */
    inline shp<BlockMatrix> getChildBT() const {return M_childBT;}

    /*! \brief Getter of the child's coupling matrix.
     *
     * \return Child's coupling matrix.
     */
    inline shp<BlockMatrix> getChildB() const {return M_childB;}

    inline void setAsInlet() {M_isInlet=true;}

    inline void setAsOutlet() {M_isOutlet=true;}

protected:
    /*! \brief Generate a quadrature rule.
     *
     * \param tag String containing the tag of the rule (e.g., STRANG10).
     * \return The quadrature rule.
     */
    shp<LifeV::QuadratureRule> generateQuadratureRule(std::string tag) const;

    /*! \brief Build the coupling vectors (columns of the coupling matrix).
     *
     * \param bfs Basis functions.
     * \param face The GeometricFace.
     * \param assembler Shared pointer to the assembler.
     *
     * \return Vector of shared pointers to the generated DistributedVectors for the coupling
     */
    std::vector<shp<DistributedVector>> buildCouplingVectors(shp<BasisFunctionFunctor> bfs,
                                                             const GeometricFace& face,
                                                             shp<aAssembler> assembler) const;

    /*! \brief Build the DOF map at the rings of the interface
     *
     * \return interfacePtr Shared pointer to a DofInterface3Dto3D object, defining the DOF map at the interface ring
     *//*
     InterfacePtrType buildRingInterfaceMap();*/

     /*! \brief Assemble the vectors for the strong coupling of velocity at the rings
      *
      * \return Vector of shared pointers to the generated DistributedVectors for the strong ring coupling
      */
     std::pair<std::vector<shp<DistributedVector>>, std::vector<shp<DistributedVector>>> buildStrongRingCouplingVectors();

    /*! \brief Identifies the interface ring DOFs employed for the strong coupling
     *
     * \return Vector of integers, corresponding to the IDs of the DOFs employed for the strong coupling
     *//*
     std::vector<LifeV::ID> identifyValidRingDOFs();*/

    /*! \brief Method to compute the coordinates of the ring mesh points in the father and child blocks
    *
    */
     void findRingPointsCoordinates();

    /*! \brief Method to identify the correspondences between father and child DOFs at the ring
    *
    */
     void buildRingDOFsMap();

    /*! \brief Evaluate in (x,y,z) a quartic RBF with center in (xc,yc,zc) and radius R
    *
    * \return Evaluation of the RBF in (x,y,z)
    */
    double evaluate_RBF(double x, double y, double z, double xc, double yc, double zc, double R);

    /// Not supported at the moment.
    std::vector<shp<DistributedVector>> buildStabilizationVectorsVelocity(shp<BasisFunctionFunctor> bfs,
                                                                          const GeometricFace& face,
                                                                          shp<aAssembler> assembler) const;

    /// Not supported at the moment.
    std::vector<shp<DistributedVector>> buildStabilizationVectorsPressure(shp<BasisFunctionFunctor> bfs,
                                                                          const GeometricFace& face,
                                                                          shp<aAssembler> assembler) const;

    /// Not supported at the moment.
    std::vector<shp<DistributedVector>> buildStabilizationVectorsLagrange() const;

    Interface                                 M_interface;

    shp<BlockMatrix>                          M_fatherBT;
    shp<BlockMatrix>                          M_fatherB;
    shp<BlockMatrix>                          M_childBT;
    shp<BlockMatrix>                          M_childB;
    shp<BlockMatrix>                          M_childBfe;
    shp<BlockMatrix>                          M_fatherBfe;

    // this is required in the RB setting to impose weakly dirichlet conditions
    shp<BlockMatrix>                          M_childBEp;
    shp<BlockMatrix>                          M_stabChild;
    shp<BlockMatrix>                          M_stabFather;

    DataContainer                             M_data;
    shp<const LifeV::MapEpetra>               M_mapLagrange;
    double                                    M_stabilizationCoupling;
    bool                                      M_isInlet;
    bool                                      M_isOutlet;
    bool                                      M_addNoSlipBC;

    std::map<LifeV::ID, Vector3D>             M_fatherRingPoints;
    std::map<LifeV::ID, Vector3D>             M_childRingPoints;
    std::map<LifeV::ID, LifeV::ID>            M_ringDOFsMap;
};

}

#endif // INTERFACEASSEMBLER_HPP
