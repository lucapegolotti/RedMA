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

#ifndef BLOCKASSEMBLER_HPP
#define BLOCKASSEMBLER_HPP

#include <redma/RedMA.hpp>
#include <redma/assemblers/abstract/aAssembler.hpp>
#include <redma/assemblers/abstract/aAssemblerFE.hpp>
#include <redma/assemblers/abstract/aAssemblerRB.hpp>
#include <redma/assemblers/AssemblerFactory.hpp>
#include <redma/assemblers/coupling/InterfaceAssembler.hpp>
#include <redma/assemblers/coupling/InletInflowAssembler.hpp>

#include <redma/reduced_basis/RBBasesManager.hpp>

#include <redma/geometry/TreeStructure.hpp>

namespace RedMA
{

/*! \brief Block assembler class to handle multiple assemblers.
 *
 * Internally, the assemblers are divided into primal assemblers, which discretize
 * the PDEs in the interior of the domains, and dual assemblers, which take care
 * of the coupling.
 */
class BlockAssembler : public aAssembler
{
    typedef DefaultAssemblersLibrary  DefaultAssemblers;

public:
    /// Default empty constructor.
    BlockAssembler() {}

    /*! \brief Constructor with multiple arguments.
     *
     * \param data DataContainer of the problem.
     * \param tree Geometric tree.
     * \param Default assemblers for the building blocks in the reference configuration.
     */
    BlockAssembler(const DataContainer& data, const TreeStructure& tree,
                   shp<DefaultAssemblers> defAssemblers = nullptr);

    /*! \brief Setup method.
     *
     * The setup function is called on the primal assemblers. Moreover,
     * the dual assemblers on the interfaces are generated based on the domain
     * connectivity within M_tree. If the RB method is used, the bases are loaded.
     */
    virtual void setup() override;

    /*! \brief Export the solutions in all subdomains.
     *
     * The export function is called in all the primal assemblers.
     *
     * \param time Current time.
     * \param sol Current solution.
     */
    virtual void exportSolution(const double& time,
                                const shp<aVector>& sol) override;

    /*! \brief Postprocess function to be called at the end of the timestep.
     *
     * The postProcess function is called in all the primal assemblers.
     *
     * \param time Current time.
     * \param sol Current solution.
     */
    virtual void postProcess(const double& time,
                             const shp<aVector>& sol) override;

    /*! \brief Getter for the mass matrix.
    *
    * The getMass method is called in all the subdomains and each output is
    * positioned in the diagonal.
    *
    * \param time Current time.
    * \param sol Current solution.
    * \return Shared pointer to BlockMatrix containing the global mass matrix.
    */
    virtual shp<aMatrix> getMass(const double& time,
                                 const shp<aVector>& sol) override;

    /*! \brief Getter for the mass matrix jacobian.
    *
    * The getMassJacobians method is called in all the subdomains and each output is
    * positioned in the diagonal.
    *
    * \param time Current time.
    * \param sol Current solution.
    * \return Shared pointer to BlockMatrix containing the global mass matrix jacobian.
    */
    virtual shp<aMatrix> getMassJacobian(const double& time,
                                         const shp<aVector>& sol) override;

    /*! \brief Getter for the global right-hand side.
    *
    * The getRightHandSide method is called in all the subdomains.
    *
    * \param time Current time.
    * \param sol Current solution.
    * \return Shared pointer to BlockVector containing the global right-hand side.
    */
    virtual shp<aVector> getRightHandSide(const double& time,
                                          const shp<aVector>& sol) override;

     /*! \brief Getter for the global right-hand side jacobian.
     *
     * The getJacobianRightHandSide method is called in all the subdomains, and the
     * coupling matrices are positioned according to the connectivity of the
     * primal assemblers.
     *
     * \param time Current time.
     * \param sol Current solution.
     * \return Shared pointer to BlockVector containing the global right-hand side.
     */
    virtual shp<aMatrix> getJacobianRightHandSide(const double& time,
                                                  const shp<aVector>& sol) override;

    /*! \brief Getter for the lifting at a specific time.
     *
     * The getLifting function is called in every primal assembler.
     * \param time Current time.
     * \return Shared pointer to BlockVector containing the lifting.
     */
    virtual shp<aVector> getLifting(const double& time) const override;

    /*! \brief Getter for the zero vector.
     *
     * The getZeroVector is called in every primal assembler.
     * \return Shared pointer to BlockVector containing the lifting.
     */
    virtual shp<aVector> getZeroVector() const override;

    /*! \brief Apply Dirichlet bcs to a matrix with prescribed
     *         diagonal coefficient.
     *
     * The applyDirichletBCsMatrix function is called in every subdomain.
     * Note: by construction, the coupling matrices have the bcs
     * incorporated in them (i.e., the rows corresponding to Dirichlet nodes are
     * zero).
     *
     * \param matrix The matrix to modify.
     * \param diagCoeff diagonal coefficient.
     */
    virtual void applyDirichletBCsMatrix(shp<aMatrix> matrix, double diagCoeff) const override;

    /*! \brief Apply homogeneous Dirichlet bcs to a vector.
     *
     * The apply0DirichletBCs function is called in every subdomain.
     *
     * \param vector The vector to modify.
     */
    virtual void apply0DirichletBCs(shp<aVector> vector) const override;

    /*! \brief Apply Dirichlet bcs to a vector at a specific time.
     *
     * The applyDirichletBCs function is called in every subdomain.
     *
     * \param time Current time
     * \param vector Shared pointer to the vector to modify.
     */
    virtual void applyDirichletBCs(const double& time, shp<aVector> vector) const override;

    /// Set exporter (set exporters in all the subdomains).
    virtual void setExporter() override;

    /*! \brief Check the magnitude of the stabilization term for the coupling.
     *
     * Currently not implemented.
     */
    virtual void checkStabTerm(const shp<aVector>& sol) const;

    /// Create and get a map of type [key,value] = [id,meshname].
    std::map<unsigned int, std::string> getIDMeshTypeMap() const;

    /*! \brief Getter for a primal assembler with given index.
     *
     * \param index Index of the primal assembler to get.
     * \return Shared pointer to aAssembler of the desired primal assembler.
     */
    inline shp<aAssembler> block(const unsigned int& index) {return M_primalAssemblers[index];}

    /*! \brief Get the primal assembler map.
     *
     * \return Retrun a map with key == id; value is a shared pointer to an
     *         aAssembler of a primal block.
     */
    std::map<unsigned int, shp<aAssembler>> getAssemblersMap() const {return M_primalAssemblers;}

    /*! \brief Get vector of dual assemblers.
     *
     * \return Vector of shared pointers to InterfaceAssembler.
     */
    std::vector<shp<InterfaceAssembler>> getDualAssemblers() const {return M_dualAssemblers;}

    /*! \brief Convert RB function to FE function.
     *
     * The convertFunctionRBtoFEM function is called on every primal assembler.
     *
     * \param rbFunction Shared pointer to BlockVector containing the RB coefficients.
     * \param comm MPI Communicator.
     * \return Shared pointer to the converted vector.
     */
    shp<aVector> convertFunctionRBtoFEM(shp<aVector> rbFunction, EPETRACOMM comm) const;

    /*! \brief Set extrapolated solution.
     *
     * The extrapolated solution is set in every primal block.
     *
     * \param exSol Shared pointer to a BlockVector containing the extrapolated
     *              solution.
     */
    virtual void setExtrapolatedSolution(const shp<aVector>& exSol) override;

    /*! \brief Getter for the nonlinear term.
     *
     * The getNonLinearTerm is called in every primal assembler.
     *
     * \return Shared pointer to a BlockVector containing the nonlinear term.
     */
    virtual shp<aVector> getNonLinearTerm() override;

    /*! \brief Get map with all the randomizible geometrical parameters.
     *
     * The map is generated by calling the same method on the TreeNodes internal
     * to the primal assemblers.
     * \return Map with key == id and value == vector of randomizible parameters.
     */
    std::map<unsigned int,std::vector<double>> getRandomizibleParametersVectors();

    /*! \brief Applies the piola transformation (or its inverse) to a BlockVector.
     *
     * \param solution Shared pointer to GlobalVector to transform.
     * \param inverse If true, inverse of Piola transformation is applied.
     */
    virtual void applyPiola(shp<aVector> solution, bool inverse) override;

    /// Initialize all the finite element spaces in the primal assemblers.
    virtual void initializeFEspaces() override;

    /*! Setter fot the default assemblers.
     *
     * \param Shared pointer to the default assemblers.
     */
    void setDefaultAssemblers(shp<DefaultAssemblers> defAssemblers) override;

    /*! Check if all the primal assemblers are finite element assemblers.
     *
     * \return True if all the primal assemblers are finite element assemblers.
     */
    bool arePrimalAssemblersFE();

    /*! Check if all the primal assemblers are reduced basis assemblers.
     *
     * \return True if all the primal assemblers are reduced basis assemblers.
     */
    bool arePrimalAssemblersRB();

protected:
    GetPot                                                        M_datafile;
    TreeStructure                                                 M_tree;
    std::map<unsigned int, shp<aAssembler>>                       M_primalAssemblers;
    std::vector<shp<InterfaceAssembler>>                          M_dualAssemblers;
    unsigned int                                                  M_numberBlocks;
    // shp<MDEIMManager>                                             M_mdeimManager;
    shp<RBBasesManager>                                           M_basesManager;
};

}

#endif // BLOCKASSEMBLER_HPP
