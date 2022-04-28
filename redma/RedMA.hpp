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

#ifndef __has_include
  static_assert(false, "__has_include not supported");
#else
#if __has_include(<filesystem>)
    #include <filesystem>
        namespace fs = std::filesystem;
    #elif __has_include(<experimental/filesystem>)
    #include <experimental/filesystem>
        namespace fs = std::experimental::filesystem;
    #endif
#endif

#include <redma/utils/PrintLog.hpp>
#include <redma/utils/Chrono.hpp>
#include <redma/utils/Exception.hpp>

#include <lifev/core/fem/FESpace.hpp>
#include <lifev/eta/fem/ETFESpace.hpp>
#include <lifev/core/array/MapEpetra.hpp>
#include <lifev/core/mesh/RegionMesh.hpp>
#include <lifev/core/util/LifeChrono.hpp>
#include <lifev/core/array/MatrixEpetra.hpp>
#include <lifev/core/array/VectorEpetra.hpp>

#include <Epetra_SerialDenseSolver.h>
#include <Epetra_SerialDenseMatrix.h>
#include <Epetra_SerialDenseVector.h>

// we define the namespace
namespace RedMA
{

class DistributedVector;
class SparseMatrix;
class DenseVector;
class DenseMatrix;

template <class Type>
using shp = std::shared_ptr<Type>;

template <class Type1, class Type2>
shp<Type1> spcast(shp<Type2> ptr)
{
    return std::static_pointer_cast<Type1>(ptr);
}

template <class Type1, class Type2>
shp<Type1> dpcast(shp<Type2> ptr)
{
    return std::dynamic_pointer_cast<Type1>(ptr);
}

}

#define COMMA                 ,

#define EPETRACOMM          shp<Epetra_Comm>

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
#define VECTOREPETRA        LifeV::VectorEpetra
#define DENSEVECTOR         Epetra_SerialDenseVector
#define MATRIXEPETRA        LifeV::MatrixEpetra<double>
#define DENSEMATRIX         Epetra_SerialDenseMatrix

#endif  // REDMA_HPP
