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

#ifndef VMSSUPGSTABILIZATION_HPP
#define VMSSUPGSTABILIZATION_HPP

#include <AbstractAssembler.hpp>
#include <Exception.hpp>

#include <lifev/eta/fem/ETFESpace.hpp>
#include <lifev/core/fem/BCHandler.hpp>
#include <lifev/eta/utils/Functions.hpp>

#include <Exception.hpp>

#include <functional>

#include <lifev/core/filter/Exporter.hpp>
#include <lifev/core/filter/ExporterVTK.hpp>
#include <lifev/core/filter/ExporterHDF5.hpp>
#include <boost/filesystem.hpp>

namespace RedMA
{

class VMS_SUPGStabilization
{
protected:
    typedef LifeV::RegionMesh<LifeV::LinearTetra>       Mesh;
    typedef std::shared_ptr<Mesh>                       MeshPtr;
    typedef std::shared_ptr<Epetra_Comm>                commPtr_Type;
    typedef LifeV::MapEpetra                            MapEpetra;
    typedef std::shared_ptr<MapEpetra>                  MapEpetraPtr;
    typedef LifeV::VectorSmall<3>                       Vector3D;
    typedef LifeV::MatrixSmall<3, 3>                    Matrix3D;
    typedef LifeV::FESpace<Mesh, MapEpetra>             FESpace;
    typedef std::shared_ptr<FESpace>                    FESpacePtr;
    typedef LifeV::ETFESpace<Mesh, MapEpetra, 3, 3>     ETFESpaceVelocity;
    typedef std::shared_ptr<ETFESpaceVelocity>          ETFESpaceVelocityPtr;
    typedef LifeV::ETFESpace<Mesh, MapEpetra, 3, 1>     ETFESpacePressure;
    typedef std::shared_ptr<ETFESpacePressure>          ETFESpacePressurePtr;
    typedef LifeV::VectorEpetra                         Vector;
    typedef std::shared_ptr<Vector>                     VectorPtr;
    typedef LifeV::MatrixEpetra<double>                 Matrix;
    typedef std::shared_ptr<Matrix>                     MatrixPtr;
public:
    VMS_SUPGStabilization(unsigned int        order,
                         FESpacePtr           fespaceVelocity,
                         FESpacePtr           fespacePressure,
                         ETFESpaceVelocityPtr etfespaceVelocity,
                         ETFESpacePressurePtr etfespacePressure);

    void setDensityAndViscosity(double density, double viscosity);

    void assembleBlocks(VectorPtr velocity, VectorPtr pressure,
                        VectorPtr velocityRhs, double dt);

private:
    double                          M_density;
    double                          M_viscosity;
    unsigned int                    M_order;
    FESpacePtr                      M_velocityFESpace;
    FESpacePtr                      M_pressureFESpace;
    ETFESpaceVelocityPtr            M_velocityFESpaceETA;
    ETFESpacePressurePtr            M_pressureFESpaceETA;
    MatrixPtr                       M_block00;
    MatrixPtr                       M_block10;
    MatrixPtr                       M_block01;
    MatrixPtr                       M_block11;
    MatrixPtr                       M_blockMass00;
    MatrixPtr                       M_blockMass10;
    MatrixPtr                       M_blockMass01;
    MatrixPtr                       M_blockMass11;
    double                          M_C_I;
};

}  // namespace RedMA

#endif  // VMSSUPGSTABILIZATION_HPP
