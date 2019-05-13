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

#ifndef NAVIERSTOKESASSEMBLER_HPP
#define NAVIERSTOKESASSEMBLER_HPP

#include <AbstractAssembler.hpp>
#include <Exception.hpp>

#include <lifev/eta/fem/ETFESpace.hpp>
#include <lifev/core/fem/BCHandler.hpp>

#include <functional>

#include <boost/filesystem.hpp>

namespace RedMA
{

class NavierStokesAssembler : public AbstractAssembler
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
    typedef LifeV::ExporterVTK<Mesh>                    Exporter;
    typedef std::shared_ptr<Exporter>                   ExporterPtr;

public:
    NavierStokesAssembler(const GetPot& datafile, commPtr_Type comm,
                          const TreeNodePtr& treeNode, bool verbose = false);

    void setup();

    MatrixPtr getMassMatrix(const unsigned int& blockrow,
                            const unsigned int& blockcol);

    MatrixPtr getJacobian(const unsigned int& blockrow,
                          const unsigned int& blockcol);

    inline unsigned int numberOfBlocks() {return 2;}

    inline unsigned int numberOfComponents() {return 3;}

    void setTimeAndPrevSolution(const double& time,
                                std::vector<VectorPtr> solution);

    void setMaxVelocityLawInflow(std::function<double(double)> maxLaw);

    void setMaxVelocityDtLawInflow(std::function<double(double)> maxLawDt);

    std::vector<VectorPtr> computeF();

    std::vector<VectorPtr> computeFder();

    void applyBCsRhsRosenbrock(std::vector<VectorPtr> rhs,
                               std::vector<VectorPtr> utilde,
                               const double& time,
                               const double& dt,
                               const double& alphai,
                               const double& gammai);

    void applyBCsBackwardEuler(std::vector<VectorPtr> rhs, const double& coeff,
                               const double& time);

    void applyBCsMatrix(MatrixPtr matrix, const double& diagonalCoefficient,
                        const unsigned int& iblock, const unsigned int& jblock);

    void setExporter();

    void exportSolutions(const double& time, std::vector<VectorPtr> solutions);

protected:
    void assembleConstantMatrices();

    void assembleStiffnessMatrix();

    void assembleDivergenceMatrix();

    void assembleMassMatrix();

    void assembleConvectiveMatrix();

    void assembleJacobianConvectiveMatrix();

    void assembleForcingTerm();

    void assembleForcingTermTimeDerivative();

    BoundaryConditionPtr createBCHandler(std::function<double(double)> law);

    static double fZero(const double& t,
                        const double& x,
                        const double& y,
                        const double& z,
                        const unsigned int& i);

    static double poiseulleInflow(const double& t,
                                  const double& x,
                                  const double& y,
                                  const double& z,
                                  const unsigned int& i,
                                  const GeometricFace& face,
                                  std::function<double(double)> maxLaw);

    void updateBCs(BoundaryConditionPtr bcToUpdate, FESpacePtr fespace);

    FESpacePtr                      M_velocityFESpace;
    FESpacePtr                      M_pressureFESpace;
    ETFESpaceVelocityPtr            M_velocityFESpaceETA;
    ETFESpacePressurePtr            M_pressureFESpaceETA;

    MatrixPtr                       M_A;
    MatrixPtr                       M_B;
    MatrixPtr                       M_Bt;
    MatrixPtr                       M_M;
    MatrixPtr                       M_C;
    MatrixPtr                       M_J;
    VectorPtr                       M_forcingTerm;
    VectorPtr                       M_forcingTermTimeDer;
    std::vector<VectorPtr>          M_prevSolution;
    double                          M_time;

    FunctionType                    M_forceFunction;
    FunctionType                    M_forceTimeDerFunction;
    std::function<double(double)>   M_maxVelocityLaw;
    std::function<double(double)>   M_maxVelocityDtLaw;

    ExporterPtr                     M_exporter;
    VectorPtr                       M_velocityExporter;
    VectorPtr                       M_pressureExporter;
};

}  // namespace RedMA

#endif  // NAVIERSTOKESASSEMBLER_HPP
