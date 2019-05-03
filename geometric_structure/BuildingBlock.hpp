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
#include <sys/stat.h>
#include <sys/types.h>
#include <iostream>

#include <Exception.hpp>
#include <PrintLog.hpp>
#include <GeometricParametersHandler.hpp>

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
#include <lifev/core/array/MatrixSmall.hpp>
#include <lifev/core/filter/ExporterHDF5.hpp>
#include <lifev/core/filter/ExporterVTK.hpp>
#include <lifev/core/fem/FESpace.hpp>

#include <boost/filesystem.hpp>

namespace RedMA
{

class GeometricFace
{
public:
    typedef LifeV::VectorSmall<3>   Vector3D;

    GeometricFace();

    GeometricFace(Vector3D center, Vector3D normal, double radius);

    void print();

    Vector3D M_center;
    Vector3D M_normal;
    double M_radius;
};

class BuildingBlock
{
protected:
    typedef LifeV::RegionMesh<LifeV::LinearTetra>          mesh_Type;
    typedef std::shared_ptr<mesh_Type>                     meshPtr_Type;
    typedef std::shared_ptr<Epetra_Comm>                   commPtr_Type;
    typedef LifeV::MapEpetra                               map_Type;
    typedef std::shared_ptr<map_Type>                      mapPtr_Type;
    typedef LifeV::VectorSmall<3>                          Vector3D;
    typedef LifeV::MatrixSmall<3,3>                        Matrix3D;
    typedef LifeV::ExporterVTK<mesh_Type>                  Exporter;
    typedef LifeV::FESpace<mesh_Type, map_Type>            FESpace_Type;
    typedef std::shared_ptr<FESpace_Type>                  FESpacePtr_Type;
    typedef LifeV::VectorEpetra                            vector_Type;
    typedef std::shared_ptr<vector_Type>                   vectorPtr_Type;
    typedef LifeV::MeshUtility::MeshTransformer<mesh_Type> Transformer;

public:
    BuildingBlock(commPtr_Type comm, bool verbose);

    void setParameterValue(std::string key, double value);

    std::map<std::string,std::shared_ptr<GeometricParameter> >&
                                                             getParametersMap();

    int readMesh(std::string meshdir = "../geometries/");

    virtual inline unsigned int expectedNumberOfChildren() = 0;

    std::string name();

    void applyAffineTransformation();

    virtual void applyNonAffineTransformation() = 0;

    void applyGlobalTransformation();

    void dumpMesh(std::string outdir, std::string meshdir,
                  std::string outputName);

    GeometricFace getOutlet(unsigned int indexFace) const;

    std::vector<GeometricFace> getOutlets() const;

    GeometricFace getInlet() const;

    void mapChildInletToParentOutlet(GeometricFace parentOutlet);

    void setIsChild(bool isChild);

    meshPtr_Type getMesh();

    void setDatafile(const GetPot& datafile);

    void setRandom();

protected:
    void applyAffineTransformationGeometricFace(GeometricFace& face,
                                                const Matrix3D& affineMatrix,
                                                const Vector3D& translation,
                                                const double& scale);

    Matrix3D computeRotationMatrix(unsigned int axis, double angle);

    Matrix3D computeRotationMatrix(Vector3D axis, double angle);

    static  void rotationFunction(double& x, double& y, double& z,
                                  const Matrix3D& affMatrix,
                                  const Vector3D& transl, const double& scale);

    static double fZero(const double& t, const double& x, const double& y,
                      const double& z, const LifeV::ID& i);

    GeometricParametersHandler M_parametersHandler;
    std::string M_name;

    std::string M_meshName;
    meshPtr_Type M_mesh;

    commPtr_Type M_comm;

    bool M_verbose;

    std::string M_datafileName;

    GeometricFace M_inlet;
    std::vector<GeometricFace> M_outlets;

    bool M_isChild;

    double M_inletScale;

    Vector3D M_inletTranslation;
    Vector3D M_inletRotationAxis;
    double M_inletAngle;

    GetPot M_datafile;
};

}  // namespace RedMA

#endif  // BUILDINGBLOCK_HPP
