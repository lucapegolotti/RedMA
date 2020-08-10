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

#ifndef DEFAULTASSEMBLERSLIBRARY_HPP
#define DEFAULTASSEMBLERSLIBRARY_HPP

#include <redma/RedMA.hpp>
#include <redma/problem/DataContainer.hpp>

#include <redma/geometry/TreeStructure.hpp>

namespace RedMA
{

class aAssembler;

SHP(aAssembler)
AssemblerFactory(const DataContainer& data, SHP(TreeNode) treeNode);

class DefaultAssemblersLibrary
{
//     typedef aAssembler   AssemblerType;
// public:
//     DefaultAssemblersLibrary(const DataContainer& data, const std::set<std::string>& meshes, EPETRACOMM comm);
//
//     SHP(TreeNode) generateDefaultTreeNode(const std::string& nameMesh);
//
//     SHP(TreeNode) generateDefaultTube(const std::string& nameMesh);
//
//     SHP(TreeNode) generateDefaultSymmetricBifurcation(const std::string& nameMesh);
//
//     SHP(AssemblerType) getDefaultAssembler(const std::string& namemesh);
//
// private:
//     DataContainer                                   M_data;
//     std::map<std::string, SHP(AssemblerType)>       M_assemblersMap;
//     EPETRACOMM                                      M_comm;
//     unsigned int                                    M_count;
};

}

#endif // DEFAULTASSEMBLERSLIBRARY_HPP
