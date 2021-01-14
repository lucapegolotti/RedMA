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

#define THRESHOLDSTAB       1e-15

namespace RedMA
{

/*! \brief Class representing an interface between two subdomains.
 *
 */
class Interface
{
    typedef aAssembler         AssemblerType;
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
};

/*! \brief Class for the assembly of coupling matrices.
 *
 */
class InterfaceAssembler
{
    typedef aAssembler         AssemblerType;

public:

    /*! \brief Constructor taking a DataContainer as argument.
     *
     * \param data The DataContainer of the problem.
     */
    InterfaceAssembler(const DataContainer& data);

    /*! \brief Constructor taking a DataContainer and an Interface as arguments.
     *
     * \param data The DataContainer of the problem.
     * \param interface The interface.
     */
    InterfaceAssembler(const DataContainer& data,
                       const Interface& interface);

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
     */
    void buildCouplingMatrices(shp<AssemblerType> assembler,
                               const GeometricFace& face,
                               shp<BlockMatrix> matrixT,
                               shp<BlockMatrix> matrix);

    /// Build stabilization matrix. Currently not implemented.
    void buildStabilizationMatrix(shp<AssemblerType> assembler,
                                  const GeometricFace& face,
                                  shp<BlockMatrix> matrix);

    /*! \brief Build Epetra map for the Lagrange multiplier space.
     *
     * \param bfs Shared pointer to the BasisFunctionFunctor.
     */
    void buildMapLagrange(shp<BasisFunctionFunctor> bfs);

    /* \brief Add coupling contribution to a right-hand side.
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

    /* \brief Add coupling contribution to the jacobian right-hand side.
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

    /* \brief Getter for the zero-vector.
     *
     * \return Shared pointer to BlockVector containing the zero vector.
     */
    shp<BlockVector> getZeroVector() const;

    /* \brief Check norm of the stabilization term.
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
     */
    std::vector<shp<DistributedVector>> buildCouplingVectors(shp<BasisFunctionFunctor> bfs,
                                                             const GeometricFace& face,
                                                             shp<aAssembler> assembler) const;

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

    Interface                              M_interface;
    shp<BlockMatrix>                       M_identity;
    shp<BlockMatrix>                       M_fatherBT;
    shp<BlockMatrix>                       M_fatherB;
    shp<BlockMatrix>                       M_childBT;
    shp<BlockMatrix>                       M_childB;
    shp<BlockMatrix>                       M_childBTfe;
    shp<BlockMatrix>                       M_childBfe;
    // this is required in the RB setting to impose weakly dirichlet conditions
    shp<BlockMatrix>                       M_childBEp;
    shp<BlockMatrix>                       M_stabChild;
    shp<BlockMatrix>                       M_stabFather;
    DataContainer                          M_data;
    shp<const LifeV::MapEpetra>            M_mapLagrange;
    double                                 M_stabilizationCoupling;
    bool                                   M_isInlet;
};

}

#endif // INTERFACEASSEMBLER_HPP
