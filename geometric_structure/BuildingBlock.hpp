// Reduced Modeling of Arteries (ReMA)
// Copyright (C) 2019  Luca Pegolotti
//
// ReMA is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// ReMA is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with this program.  If not, see <http://www.gnu.org/licenses/>.

#ifndef BUILDINGBLOCK_HPP
#define BUILDINGBLOCK_HPP

#include <map>
#include <memory>

#include <Exception.hpp>
#include <PrintLog.hpp>

#include <lifev/core/mesh/RegionMesh.hpp>
#include <lifev/core/mesh/MeshData.hpp>
#include <lifev/core/mesh/MeshPartitioner.hpp>
#include <lifev/core/filter/GetPot.hpp>

#include <Epetra_MpiComm.h>

namespace ReMA
{

class BuildingBlock
{
protected:
    typedef LifeV::RegionMesh<LifeV::LinearTetra>   mesh_Type;
    typedef std::shared_ptr<mesh_Type>              meshPtr_Type;

    typedef std::shared_ptr<Epetra_Comm>            commPtr_Type;

public:
    BuildingBlock(commPtr_Type comm, bool verbose);

    void setParameterValue(std::string key, double value);

    int readMesh(std::string meshdir);

protected:
    std::map<std::string,double> M_parametersMap;
    std::string M_name;

    std::string M_meshName;
    meshPtr_Type M_mesh;

    commPtr_Type M_comm;

    bool M_verbose;

    std::string M_datafileName;
};

}  // namespace ReMA

#endif  // BUILDINGBLOCK_HPP
