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

#ifndef GLOBALSOLVER_HPP
#define GLOBALSOLVER_HPP

#include <GlobalAssembler.hpp>
#include <TreeStructure.hpp>
#include <GeometryParser.hpp>
#include <TreeStructure.hpp>

#include <lifev/core/array/MapVector.hpp>
#include <lifev/core/array/MapEpetra.hpp>

namespace RedMA
{

template <class AssemblerType>
class GlobalSolver
{
    typedef LifeV::MatrixEpetraStructured<double>           MatrixStructured;
    typedef std::shared_ptr<Epetra_Comm>                    commPtr_Type;
    typedef LifeV::MapEpetra                                map_Type;
	typedef std::shared_ptr<map_Type>                       mapPtr_Type;
    typedef LifeV::MapVector<map_Type>                      MapVector;
    typedef std::shared_ptr<MapVector>                      MapVectorPtr;
    typedef std::shared_ptr<MatrixStructured>               MatrixStructuredPtr;

public:
    GlobalSolver(std::string xmlFile, std::string geometriesDir,
                 commPtr_Type comm, bool verbose);
private:
    GeometryParser M_geometryParser;
    TreeStructure M_tree;
    MapVectorPtr M_mapVector;
    MatrixStructuredPtr M_globalMatrix;
};

}  // namespace RedMA

#include "GlobalSolver_imp.hpp"

#endif  // GLOBALSOLVER_HPP
