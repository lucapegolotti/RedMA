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
#include <Teuchos_ParameterList.hpp>
#include <Teuchos_XMLParameterListHelpers.hpp>
#include <Teuchos_RCP.hpp>

namespace RedMA
{

class NonAffineDeformer
{
    typedef LifeV::RegionMesh<LifeV::LinearTetra>          mesh_Type;
    typedef shp<mesh_Type>                     meshPtr_Type;
    typedef shp<Epetra_Comm>                   commPtr_Type;
    typedef LifeV::MapEpetra                               map_Type;
    typedef shp<map_Type>                      mapPtr_Type;
    typedef LifeV::FESpace<mesh_Type, map_Type>            fespace_Type;
    typedef shp<fespace_Type>                  fespacePtr_Type;
    typedef LifeV::ETFESpace< mesh_Type, map_Type, 3, 3 >  fespaceETA_Type;
    typedef shp<fespaceETA_Type>               fespaceETAPtr_Type;
    typedef LifeV::MatrixEpetra<LifeV::Real>               matrix_Type;
    typedef shp<matrix_Type>                   matrixPtr_Type;
    typedef shp<LifeV::BCHandler>              bcPtr_Type;
    typedef LifeV::VectorEpetra                            vector_Type;
    typedef shp<vector_Type>                   vectorPtr_Type;


public:
    NonAffineDeformer(meshPtr_Type mesh, commPtr_Type comm,
                      bool verbose = false);

    void setXMLsolver(std::string filename);

    void applyBCs(bcPtr_Type bcs);

    void deformMesh(LifeV::MeshUtility::MeshTransformer<mesh_Type>& transformer);

private:
    void assembleStiffness(const double young, const double poisson);

    vectorPtr_Type solveSystem();

    meshPtr_Type M_mesh;
    commPtr_Type M_comm;
    bool M_verbose;

    fespacePtr_Type M_fespace;
    fespaceETAPtr_Type M_fespaceETA;

    matrixPtr_Type M_stiffness;
    vectorPtr_Type M_rhs;

    std::string M_XMLsolver;
};

}  // namespace RedMA

#endif  // NONAFFINEDEFORMER_HPP
