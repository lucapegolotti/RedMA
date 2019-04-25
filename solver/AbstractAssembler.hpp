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

#ifndef ABSTRACTASSEMBLER_HPP
#define ABSTRACTASSEMBLER_HPP

#include <TreeStructure.hpp>

#include <lifev/core/array/MapEpetra.hpp>

namespace RedMA
{

class AbstractAssembler
{
protected:
    typedef std::shared_ptr<TreeNode>                       TreeNodePtr;
    typedef LifeV::MapEpetra                                MapEpetra;
    typedef std::shared_ptr<MapEpetra>                      MapEpetraPtr;
    typedef LifeV::MapEpetra                                map_Type;
    typedef std::shared_ptr<map_Type>                       mapPtr_Type;
    typedef LifeV::MapVector<map_Type>                      MapVector;
    typedef std::shared_ptr<MapVector>                      MapVectorPtr;
    typedef std::vector<MapEpetraPtr>                       MapVectorSTD;

public:
    AbstractAssembler();

    AbstractAssembler(const TreeNodePtr& treeNode);

    void addMapsToVector(MapVectorPtr& mapVector);

protected:
    TreeNodePtr               M_treeNode;
    std::vector<MapEpetraPtr> M_maps;
};

}  // namespace RedMA

#endif  // ABSTRACTASSEMBLER_HPP
