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

#include <redma/solver/array/MatrixEp.hpp>
#include <redma/solver/array/VectorEp.hpp>

// we define the namespace
namespace RedMA
{

}

// this is a trick to allow for commas in arguments to SHP macro
template<typename T> struct argument_type;
template<typename T, typename U> struct argument_type<T(U)> {typedef U type;};

#define SHP(TYPE)           std::shared_ptr<typename argument_type<void(TYPE)>::type>

#define EPETRACOMM          SHP(Epetra_Comm)

#define FEVECTOR            RedMA::VectorEp
#define FEMATRIX            RedMA::MatrixEp

#endif  // REDMA_HPP
