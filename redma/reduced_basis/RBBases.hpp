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

#ifndef RBBASES_HPP
#define RBBASES_HPP

#include <redma/RedMA.hpp>
#include <redma/problem/DataContainer.hpp>

#include <redma/array/BlockMatrix.hpp>

#include <redma/RedMA.hpp>

namespace RedMA
{

/*! \brief Class handling a reduced basis.
 */
class RBBases
{
    typedef boost::numeric::ublas::matrix<std::vector<shp<VECTOREPETRA>>> GridVectors;

public:
    /*! \brief Constructor.
     *
     * \param data A DataContainer.
     * \param comm The MPI Communicator.
     */
    RBBases(const DataContainer& data,
            EPETRACOMM comm);

    /*! \brief Set the number of fields.
     *
     * \param nfields A DataContainer.
     * \param comm The MPI Communicator.
     */
    void setNumberOfFields(const unsigned int& nfields);

    /*! \brief Set path.
     *
     * This is the path were the reduced basis must be saved or are stored.
     *
     * \param path The path.
     */
    void setPath(std::string path);

    /// Load the singular values from file.
    void loadSingularValues();

    /// Load the singular values from file.
    void loadBases();

    /// Load the singular values from file.
    void dump();

    /*! \brief Set the finite element space of a field.
     *
     * \param fespace The finite element space.
     * \param indexbasis Index of the component.
     */
    void setFESpace(shp<FESPACE> fespace,
                    const unsigned int& indexbasis);

    /*! \brief Get basis vectors corresponding to an index (w/o supremizers).
     *
     * The number of returned vectors depends on the POD tolerance.
     *
     * \param index The index of the field.
     * \return A std vector containing the basis vectors.
     */
    std::vector<shp<VECTOREPETRA>> getBasis(const unsigned int& index);


    /*! \brief Get all basis vectors corresponding to an index (w/o supremizers).
     *
     * \param index The index of the field.
     * \return A std vector containing the basis vectors.
     */
    std::vector<shp<VECTOREPETRA>> getFullBasis(const unsigned int& index)
    {
        return M_bases[index];
    }

    /*! \brief Get basis vectors corresponding to an index (with supremizers).
     *
     * The vectors may also be scaled with the Piola transformation. Hence,
     * we need to pass the ID of the desired block.
     *
     * \param index The index of the field.
     * \param ID The ID of the block.
     * \return A std vector containing the basis vectors.
     */
    std::vector<shp<VECTOREPETRA>> getEnrichedBasis(const unsigned int& index,
                                                    unsigned int ID);

    /*! \brief Get primal supremizers of one field with respect to another.
     *
     * \param i The index of the field which must be augmented.
     * \param j The index of the field which poses the constraint.
     */
    std::vector<shp<VECTOREPETRA>> getPrimalSupremizers(const unsigned int& i,
                                                        const unsigned int& j)
    {
        return M_primalSupremizers(i,j);
    }

    /*! \brief Project a matrix on both sides onto a reduced basis.
     *
     * \param matrix The matrix to project.
     * \param basisIndexRow The index of the basis on the left.
     * \param basisIndexCol The index of the basis on the right.
     * \param ID The ID of the building block.
     * \return The projected matrix.
     */
    shp<aMatrix> matrixProject(shp<aMatrix> matrix,
                               unsigned int basisIndexRow,
                               unsigned int basisIndexCol,
                               unsigned int ID,
                               shp<aMatrix> normMatrix=nullptr);

    /*! \brief Project a block matrix on the left onto a reduced basis.
     *
     * \param matrix The matrix to project.
     * \param ID The ID of the building block.
     * \return The projected matrix.
     */
    shp<BlockMatrix> leftProject(shp<BlockMatrix> matrix,
                                 unsigned int ID);

    /*! \brief Project a sparse matrix on the left.
     *
     * \param matrix The matrix to project.
     * \param basisIndex The index of the basis.
     * \param ID The ID of the building block.
     * \return The projected matrix.
     */
    shp<DenseMatrix> leftProject(shp<SparseMatrix> matrix,
                                 unsigned int basisIndex,
                                 unsigned int ID);

    /*! \brief Project a vector on the left.
     *
     * \param vector The vector to project.
     * \param ID The ID of the building block.
     * \return The projected vector.
     */
    shp<aVector> leftProject(shp<aVector> vector,
                             unsigned int basisIndex,
                             unsigned int ID,
                             shp<aMatrix> normMatrix=nullptr);

    /*! \brief Project a block vector on the left.
     *
     * \param vector The vector to project.
     * \param ID The ID of the building block.
     * \return The projected vector.
     */
    shp<BlockVector> leftProject(shp<BlockVector> vector,
                                 unsigned int ID);

    /*! \brief Project a distributed vector on the left.
     *
     * \param vector The vector to project.
     * \param basisIndex The index of the basis
     * \param ID The ID of the building block.
     * \return The projected vector.
     */
    shp<DenseVector> leftProject(shp<DistributedVector> vector,
                                 unsigned int basisIndex,
                                 unsigned int ID);

    /*! \brief Project a vector onto the basis for the Lagrange multipliers.
     *
     * The basis of the Lagrange multiplier is the identity. Hence,
     * this function simply performs the conversion from distributed to dense.
     *
     * \param vector The vector to project.
     * \return The projected vector.
     */
    shp<BlockVector> projectOnLagrangeSpace(shp<BlockVector> vector);

    /*! \brief Project a block matrix on the right onto a reduced basis.
     *
     * \param matrix The matrix to project.
     * \param ID The ID of the building block.
     * \return The projected matrix.
     */
    shp<BlockMatrix> rightProject(shp<BlockMatrix> matrix,
                                  unsigned int ID);

    /*! \brief Project a sparse matrix on the right.
     *
     * \param matrix The matrix to project.
     * \param basisIndex The index of the basis.
     * \param ID The ID of the building block.
     * \return The projected matrix.
     */
    shp<DenseMatrix> rightProject(shp<SparseMatrix> matrix,
                                  unsigned int basisIndex,
                                  unsigned int ID);

    /*! \brief Recover finite element function from a reduced one.
     *
     * \param rbSolution The reduced solution.
     * \param index The index of the basis.
     * \param ID The ID of the building block.
     * \return The finite element solution.
     */
    shp<VECTOREPETRA> reconstructFEFunction(shp<aVector> rbSolution,
                                            unsigned int index,
                                            unsigned int ID);

    /*! \brief Add a primal supremizer to a set.
     *
     * \param supremizer The supremizer.
     * \param fieldToAugment The index of the field to augment.
     * \param fieldConstraint The index of the field posing the constraint.
     */
    void addPrimalSupremizer(shp<VECTOREPETRA> supremizer,
                             const unsigned int& fieldToAugment,
                             const unsigned int& fieldConstraint);

    /*! \brief Add a dual supremizer to a set.
     *
     * \param supremizer The supremizer.
     * \param fieldToAugment The index of the field to augment.
     */
    void addDualSupremizer(shp<VECTOREPETRA> supremizer,
                           const unsigned int& fieldToAugment);

    /*! \brief Get the size of a basis.
     *
     * \param index The index of the basis.
     * \return The size of the basis.
     */
    inline unsigned int getSizeBasis(const unsigned int& index)
    {
        return getBasis(index).size();
    }

    /*! \brief Get the full size of a basis (regardless of POD tolerance).
     *
     * \param index The index of the basis.
     * \return The size of the basis.
     */
    inline unsigned int getFullSizeBasis(const unsigned int& index)
    {
        return M_bases[index].size();
    }

    /*! \brief Get the size of an enriched basis (included supremizers).
     *
     * \param index The index of the basis.
     * \return The size of the basis.
     */
    inline unsigned int getSizeEnrichedBasis(const unsigned int& index)
    {
        return getEnrichedBasis(index,0).size();
    }


    // void normalizeBasis(const unsigned int& index,
    //                     shp<MATRIXEPETRA> normMatrix);

    /*! \brief Compute number of basis functions to use online.
     *
     * The number depends on the POD tolerance associated with the basis.
     *
     * \param index The index of the basis.
     */
    void computeOnlineNumberBasisFunctions(unsigned int index);

    /*! \brief Get components to keep in a given basis (based on POD tolerance).
     *
     * \param index The index of the basis.
     * \return The vector of components.
     */
    std::vector<unsigned int> getSelectors(unsigned int index);

    /// Print details of the reduced basis.
    void print();

    /*! \brief Scale a reduced basis with the Piola transformation.
     *
     * \param index The index of the basis.
     * \param ID The index of the building block.
     * \param transform The transformation function.
     */
    void scaleBasisWithPiola(unsigned int index,
                             unsigned int ID,
                             std::function<void(shp<VECTOREPETRA>)> transform);

    /*! \brief Set basis functions associated with an index.
     *
     * \param basisFunctions The functions to set.
     * \param index The index of the basis.
     */
    void setBasisFunctions(std::vector<shp<VECTOREPETRA>> basisFunctions,
                           unsigned int index)
    {
        M_bases[index] = basisFunctions;
    }

    /// Get the number of fields
    inline unsigned int getNumFields()
    {
        return M_numFields;
    }

    /*! \brief Check if a basis has supremizers.
     *
     * \param index Index of the basis.
     * \return If true, the basis has supremizers.
     */
    inline bool hasSupremizers(const unsigned int& index)
    {
        for (unsigned int j = 0; j < M_numFields; j++)
        {
            if (M_primalSupremizers(index,j).size() > 0)
                return true;
        }

        if (M_dualSupremizers[index].size() > 0)
            return true;

        return false;
    }

    /*! \brief Get enriched basis matrices with a certain index.
     *
     * \param index Index of the basis.
     * \param ID ID of the building block.
     * \param transpose If true, the returned matrix is transposed.
     */
    shp<SparseMatrix> getEnrichedBasisMatrices(const unsigned int& index,
                                               const unsigned int& ID,
                                               bool transpose);

private:
    void addVectorsFromFile(std::string filename,
                            std::vector<shp<VECTOREPETRA>>& vectors,
                            const unsigned int& indexField,
                            int Nmax = -1);

    unsigned int                                                      M_numFields;
    DataContainer                                                     M_data;
    EPETRACOMM                                                        M_comm;
    std::vector<std::vector<shp<VECTOREPETRA>>>                       M_bases;
    std::vector<shp<SparseMatrix>>                                    M_enrichedBasesMatrices;
    std::vector<shp<SparseMatrix>>                                    M_enrichedBasesMatricesTransposed;
    // this is a grid because the row indicates the field to be augmented (velocity)
    // and the column indicates the constraining field (pressure)
    GridVectors                                                       M_primalSupremizers;
    std::vector<std::vector<shp<VECTOREPETRA>>>                       M_dualSupremizers;
    std::string                                                       M_meshName;
    std::string                                                       M_path;
    std::vector<std::vector<double>>                                  M_svs;
    std::vector<shp<FESPACE>>                                         M_fespaces;
    bool                                                              M_fespacesAreSet;
    // max of tolerance on the individual fields
    double                                                            M_onlineTol;
    // online tolerance on the individual fields
    std::vector<double>                                               M_onlineTols;
    std::vector<unsigned int>                                         M_NsOnline;
    // these are needed only when we use a different basis in every block (e.g.
    // when we rescale with piola)
    std::map<unsigned int,std::map<unsigned int,
                          std::vector<shp<VECTOREPETRA>>>>            M_enrichedBasesMap;
    std::map<unsigned int,std::map<unsigned int,shp<SparseMatrix>>>   M_enrichedBasesMatricesMap;
    std::map<unsigned int,std::map<unsigned int,shp<SparseMatrix>>>   M_enrichedBasesMatricesTransposedMap;
};

}  // namespace RedMA

#endif  // RBBASES_HPP
