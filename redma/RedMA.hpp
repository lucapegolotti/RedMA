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
#ifndef REDMA_HPP
#define REDMA_HPP

#include <Epetra_ConfigDefs.h>
#ifdef EPETRA_MPI
#include <mpi.h>
#include <Epetra_MpiComm.h>
#else
#include <Epetra_SerialComm.h>
#endif

#include <memory>

#include <redma/array/SparseMatrix.hpp>
#include <redma/array/DistributedVector.hpp>
#include <redma/array/DenseVector.hpp>
#include <redma/array/DenseMatrix.hpp>
#include <redma/utils/PrintLog.hpp>
#include <redma/utils/Chrono.hpp>

#include <lifev/core/fem/FESpace.hpp>
#include <lifev/eta/fem/ETFESpace.hpp>
#include <lifev/core/array/MapEpetra.hpp>
#include <lifev/core/mesh/RegionMesh.hpp>
#include <lifev/core/util/LifeChrono.hpp>

#include <Epetra_SerialDenseSolver.h>

// we define the namespace
namespace RedMA
{

}

// this is a trick to allow for commas in arguments to SHP macro
// template<typename T> struct argument_type;
// template<typename T, typename U> struct argument_type<T(U)> {typedef U type;};
#define COMMA                 ,

// #define SHP(TYPE)           std::shared_ptr<argument_type<void(TYPE)>::type>
#define SHP(TYPE)           std::shared_ptr<TYPE>

#define EPETRACOMM          SHP(Epetra_Comm)

#define FEVECTOR            RedMA::DistributedVector
#define FEMATRIX            RedMA::SparseMatrix
#define RBVECTOR            RedMA::DenseVector
#define RBMATRIX            RedMA::DenseMatrix
#define DENSESOLVER         Epetra_SerialDenseSolver

#define MESH                LifeV::RegionMesh<LifeV::LinearTetra>
#define MAPEPETRA           LifeV::MapEpetra
#define FESPACE             LifeV::FESpace<MESH COMMA MAPEPETRA>
#define ETFESPACE3          LifeV::ETFESpace<MESH COMMA MAPEPETRA COMMA 3 COMMA 3>
#define ETFESPACE1          LifeV::ETFESpace<MESH COMMA MAPEPETRA COMMA 3 COMMA 1>

#endif  // REDMA_HPP
