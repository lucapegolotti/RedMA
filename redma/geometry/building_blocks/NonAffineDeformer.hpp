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

#ifndef NONAFFINEDEFORMER_HPP
#define NONAFFINEDEFORMER_HPP

#include <redma/RedMA.hpp>

#include <lifev/core/array/MatrixEpetra.hpp>
#include <lifev/core/array/VectorEpetra.hpp>

#include <lifev/core/fem/BCManage.hpp>
#include <lifev/eta/fem/ETFESpace.hpp>
#include <lifev/core/fem/FESpace.hpp>
#include <lifev/eta/expression/Integrate.hpp>
#include <lifev/core/mesh/MeshUtility.hpp>
#include <lifev/core/algorithm/LinearSolver.hpp>
#include <lifev/core/algorithm/PreconditionerML.hpp>
#include <lifev/core/algorithm/PreconditionerIfpack.hpp>
#include <Teuchos_ParameterList.hpp>
#include <Teuchos_XMLParameterListHelpers.hpp>
#include <Teuchos_RCP.hpp>

namespace RedMA
{

/*! \brief Class to deform a mesh by solving a linear elasticity problem.
 *
 * The boundary conditions are imposed such that specific portions of the
 * boundary are moved in a desired position. In these regions, we impose
 * Dirichlet boundary conditions.
 */
class NonAffineDeformer
{
public:

    /*! \brief Constructor.
     *
     * \param mesh The mesh to transform.
     * \param comm The MPI Communicator.
     * \param verbose If true, output is pushed to standard output.
     */
    NonAffineDeformer(shp<MESH> mesh,
                      EPETRACOMM comm,
                      bool verbose = false);

    /*! \brief Set the file with the linear solver data.
     *
     * \param filename The name of the file.
     */
    void setXMLsolver(std::string filename);

    /*! \brief Apply boundary conditions to the linear system.
     *
     * \param bcs The boundary conditions.
     */
    void applyBCs(shp<LifeV::BCHandler> bcs);

    /*! \brief Deform the mesh.
     *
     * \param transformer The mesh transformer.
     */
    void deformMesh(LifeV::MeshUtility::MeshTransformer<MESH>& transformer);

    shp<VECTOREPETRA> solveSystem(const std::string& precType = "ML");

    void deformMeshComposite(LifeV::MeshUtility::MeshTransformer<MESH>& transformer, shp<VECTOREPETRA> displacement);

private:
    void assembleStiffness(const double young, const double poisson);

    shp<MESH>               M_mesh;
    EPETRACOMM              M_comm;
    bool                    M_verbose;

    shp<FESPACE>            M_fespace;
    shp<ETFESPACE3>         M_fespaceETA;

    shp<MATRIXEPETRA>       M_stiffness;
    shp<VECTOREPETRA>       M_rhs;

    std::string             M_XMLsolver;
};

}  // namespace RedMA

#endif  // NONAFFINEDEFORMER_HPP
