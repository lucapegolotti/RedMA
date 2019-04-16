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

#include <lifev/core/LifeV.hpp>
#include <lifev/core/array/VectorEpetra.hpp>
#include <Epetra_SerialComm.h>
#include <Epetra_MpiComm.h>
#include <lifev/core/mesh/MeshUtility.hpp>
#include <lifev/core/array/VectorSmall.hpp>
#include <lifev/core/filter/ExporterHDF5.hpp>
#include <lifev/core/filter/ExporterVTK.hpp>
#include <lifev/core/fem/FESpace.hpp>

namespace RedMA
{

class BuildingBlock
{
protected:
    typedef LifeV::RegionMesh<LifeV::LinearTetra>   mesh_Type;
    typedef std::shared_ptr<mesh_Type>              meshPtr_Type;
    typedef std::shared_ptr<Epetra_Comm>            commPtr_Type;
    typedef LifeV::MapEpetra                        map_Type;
    typedef std::shared_ptr<map_Type>               mapPtr_Type;
    typedef LifeV::VectorSmall<3>                   Vector3D;
    typedef LifeV::ExporterVTK<mesh_Type>           Exporter;
    typedef LifeV::FESpace<mesh_Type, map_Type>     FESpace_Type;
    typedef std::shared_ptr<FESpace_Type>           FESpacePtr_Type;
    typedef LifeV::VectorEpetra                     vector_Type;
    typedef std::shared_ptr<vector_Type>            vectorPtr_Type;

public:
    BuildingBlock(commPtr_Type comm, bool verbose);

    void setParameterValue(std::string key, double value);

    std::map<std::string,double>& getParametersMap();

    int readMesh(std::string meshdir);

    virtual inline unsigned int expectedNumberOfChildren() = 0;

    std::string name();

    virtual void applyAffineTransformation();

    void dumpMesh(std::string outdir, std::string meshdir,
                  std::string outputName);

protected:
    std::map<std::string,double> M_parametersMap;
    std::string M_name;

    std::string M_meshName;
    meshPtr_Type M_mesh;

    commPtr_Type M_comm;

    bool M_verbose;

    std::string M_datafileName;
};

}  // namespace RedMA

#endif  // BUILDINGBLOCK_HPP
