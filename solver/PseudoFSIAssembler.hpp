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

#ifndef PSEUDOFSIASSEMBLER_HPP
#define PSEUDOFSIASSEMBLER_HPP

#include <NavierStokesAssembler.hpp>
#include <lifev/core/fem/TimeAndExtrapolationHandler.hpp>

namespace RedMA
{

class PseudoFSIAssembler : public NavierStokesAssembler
{
protected:
    typedef LifeV::ETFESpace<Mesh, MapEpetra, 3, 3>     ETFESpaceVelocity;
    typedef std::shared_ptr<ETFESpaceVelocity>          ETFESpaceVelocityPtr;
    typedef LifeV::ETFESpace<Mesh, MapEpetra, 3, 1>     ETFESpacePressure;
    typedef std::shared_ptr<ETFESpacePressure>          ETFESpacePressurePtr;
    typedef std::function<double(double const&,
                                 double const&,
                                 double const&,
                                 double const&,
                                 unsigned int const& )> FunctionType;

    typedef std::shared_ptr<LifeV::BCHandler>           BoundaryConditionPtr;
    typedef LifeV::ExporterVTK<Mesh>                    ExporterVTK;
    typedef LifeV::ExporterHDF5<Mesh>                   ExporterHDF5;
    typedef LifeV::Exporter<Mesh>                       Exporter;
    typedef std::shared_ptr<Exporter>                   ExporterPtr;

public:
    PseudoFSIAssembler(const GetPot& datafile, commPtr_Type comm,
              const TreeNodePtr& treeNode, bool verbose = false);

    std::vector<VectorPtr> computeF() override;

    std::vector<VectorPtr> computeFder() override;

    void postProcess() override;

    void setExporter() override;

    void setup() override;

    virtual unsigned int numberOfBlocks() override {return 3;};

    void getMassMatix();

    virtual MatrixPtr getMassMatrix(const unsigned int& blockrow,
                                    const unsigned int& blockcol) override;

    virtual MatrixPtr getJacobian(const unsigned int& blockrow,
                                  const unsigned int& blockcol) override;

protected:
    void assembleMassMatrix() override;

    void assembleStiffnessMatrix() override;

    void computeLameConstants();

    void computeBoundaryIndicator();

    double                             M_lameI;
    double                             M_lameII;
    MatrixPtr                          M_boundaryStiffness;
    MatrixPtr                          M_boundaryMass;
    MatrixPtr                          M_massDisplacement;
    VectorPtr                          M_displacementExporter;
    VectorPtr                          M_boundaryIndicator;
};

}  // namespace RedMA

#endif  // PSEUDOFSIASSEMBLER_HPP
