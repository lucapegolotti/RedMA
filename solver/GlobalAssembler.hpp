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

#ifndef GLOBALASSEMBLER_HPP
#define GLOBALASSEMBLER_HPP

#include <TreeStructure.hpp>

#include <lifev/core/array/MapVector.hpp>
#include <lifev/core/array/MatrixEpetraStructured.hpp>
#include <lifev/core/array/VectorEpetraStructured.hpp>

namespace RedMA
{

template <class AssemblerType>
class GlobalAssembler
{
    typedef std::shared_ptr<TreeNode>                       TreeNodePtr;
    typedef std::shared_ptr<AssemblerType>                  AssemblerTypePtr;
    typedef LifeV::MapEpetra                                map_Type;
	typedef std::shared_ptr<map_Type>                       mapPtr_Type;
    typedef LifeV::MapVector<map_Type>                      MapVector;
    typedef std::shared_ptr<MapVector>                      MapVectorPtr;
    typedef LifeV::VectorEpetraStructured                   Vector;
    typedef std::shared_ptr<Vector>                         VectorPtr;
    typedef LifeV::MatrixEpetraStructured<double>           Matrix;
    typedef std::shared_ptr<Matrix>                         MatrixPtr;
    typedef std::shared_ptr<Epetra_Comm>                    commPtr_Type;

public:
    GlobalAssembler(const GetPot& datafile, commPtr_Type comm,
                    bool verbose = false);

    void buildPrimalStructures(TreeStructure& tree);

    MapVectorPtr getMapVector() const;

    MatrixPtr getGlobalMass() const;

    MatrixPtr assembleJacobianF(const double& time, VectorPtr u) const;

    VectorPtr computeF(const double& time, VectorPtr u) const;

    VectorPtr computeFder(const double& time, VectorPtr u) const;

private:
    std::map<unsigned int, AssemblerTypePtr> M_assemblersMap;
    GetPot                                   M_datafile;
    commPtr_Type                             M_comm;
    bool                                     M_verbose;
    MapVectorPtr                             M_mapVector;
    MatrixPtr                                M_massMatrix;
};

}  // namespace RedMA

#include "GlobalAssembler_imp.hpp"

#endif  // GLOBALASSEMBLER_HPP
