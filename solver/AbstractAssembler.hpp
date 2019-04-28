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
#include <PrintLog.hpp>

#include <lifev/core/array/MapEpetra.hpp>
#include <lifev/core/filter/GetPot.hpp>
#include <lifev/core/array/MatrixEpetra.hpp>

namespace RedMA
{

class AbstractAssembler
{
protected:
    typedef std::shared_ptr<TreeNode>                      TreeNodePtr;
    typedef LifeV::MapEpetra                               MapEpetra;
    typedef std::shared_ptr<MapEpetra>                     MapEpetraPtr;
    typedef std::vector<MapEpetraPtr>                      MapVectorSTD;
    typedef LifeV::RegionMesh<LifeV::LinearTetra>          Mesh;
    typedef std::shared_ptr<Mesh>                          MeshPtr;
    typedef std::shared_ptr<Epetra_Comm>                   commPtr_Type;
    typedef LifeV::VectorSmall<3>                          Vector3D;
    typedef LifeV::MatrixSmall<3,3>                        Matrix3D;
    typedef LifeV::FESpace<Mesh, MapEpetra>                FESpace;
    typedef std::shared_ptr<FESpace>                       FESpacePtr;
    typedef LifeV::VectorEpetra                            Vector;
    typedef std::shared_ptr<Vector>                        VectorPtr;
    typedef LifeV::MatrixEpetra<double>                    Matrix;
    typedef std::shared_ptr<Matrix>                        MatrixPtr;

public:
    AbstractAssembler(const GetPot& datafile, commPtr_Type comm,
                      const TreeNodePtr& treeNode, bool verbose = false);

    void addMaps(MapEpetraPtr& globalMap, std::vector<unsigned int>& dimensions);

    // this method should be used to:
    // 1) create the finite element spaces
    // 2) assemble the constant matrices
    virtual void setup() = 0;

    inline virtual unsigned int numberOfBlocks() = 0;

protected:
    TreeNodePtr                 M_treeNode;
    std::vector<MapEpetraPtr>   M_maps;
    GetPot                      M_datafile;
    commPtr_Type                M_comm;
    bool                        M_verbose;
};

}  // namespace RedMA

#endif  // ABSTRACTASSEMBLER_HPP
