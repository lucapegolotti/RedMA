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

#include <BlockVector.hpp>

#include <functional>

#include <lifev/core/filter/Exporter.hpp>
#include <lifev/core/filter/ExporterVTK.hpp>
#include <lifev/core/filter/ExporterHDF5.hpp>
#include <boost/filesystem.hpp>

#include <SUPGStabilization.hpp>

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
                                 unsigned int const&)>  FunctionType;

    typedef std::function<double(double const&,
                                 double const&,
                                 double const&,
                                 double const&,
                                 unsigned int const&,
                                 LifeV::VectorSmall<3> const&)> FunctionNeumannType;

    typedef std::shared_ptr<LifeV::BCHandler>           BoundaryConditionPtr;
    typedef LifeV::ExporterVTK<Mesh>                    ExporterVTK;
    typedef LifeV::ExporterHDF5<Mesh>                   ExporterHDF5;
    typedef LifeV::Exporter<Mesh>                       Exporter;
    typedef std::shared_ptr<Exporter>                   ExporterPtr;
    typedef std::shared_ptr<AbstractMatrix>             AbstractMatrixPtr;
    typedef std::shared_ptr<AbstractVector>             AbstractVectorPtr;

public:
    NavierStokesAssembler(const GetPot& datafile, commPtr_Type comm,
                          const TreeNodePtr& treeNode, bool verbose = false);

    virtual void setup() override;

    virtual AbstractMatrixPtr assembleMassMatrix(double* diagonalCoefficient = nullptr) override;

    virtual MapEpetraPtr getGlobalMap() const override {toBeImplemented();};

    virtual AbstractMatrixPtr getMassMatrix() override;

    virtual MatrixPtr getJacobian(const unsigned int& blockrow,
                                  const unsigned int& blockcol) override;

    virtual MatrixPtr getJacobianPrec(const unsigned int& blockrow,
                                      const unsigned int& blockcol) override;

    virtual unsigned int numberOfBlocks() override {return 2;}

    virtual inline unsigned int numberOfComponents() override {return 3;}

    virtual void setTimeAndPrevSolution(const double& time,
                                        BlockVector solution,
                                        bool assembleBlocks = true) override;

    virtual void setLawInflow(std::function<double(double)> maxLaw) override;

    virtual void setLawDtInflow(std::function<double(double)> maxLawDt) override;

    virtual BlockVector computeF() override;

    virtual BlockVector computeFder() override;

    virtual void applyBCsRhsRosenbrock(BlockVector rhs,
                                       BlockVector utilde,
                                       const double& time,
                                       const double& dt,
                                       const double& alphai,
                                       const double& gammai) override;

    virtual void applyBCsBackwardEuler(BlockVector rhs, const double& coeff,
                                       const double& time) override;

    virtual void applyBCsMatrix(MatrixPtr matrix, const double& diagonalCoefficient,
                                const unsigned int& iblock, const unsigned int& jblock) override;

    virtual void setExporter();

    virtual void exportSolutions(const double& time, BlockVector solutions) override;

    virtual std::vector<double> computeNorms(BlockVector solutions) override;

    virtual std::vector<double> computeErrors(BlockVector solutions,
                                      const double& time) override;

    static std::string normFileFirstLine();

    static std::string errorsFileFirstLine();

    virtual MatrixPtr getUpdateMass(const unsigned int& blockrow,
                                    const unsigned int& blockcol) override;

    virtual MatrixPtr getUpdateMassJac(const unsigned int& blockrow,
                                       const unsigned int& blockcol) override;

    virtual MatrixPtr getUpdateMassJacVelocity(const unsigned int& blockrow,
                                              const unsigned int& blockcol) override;

    BlockVector initialCondition() override;

    virtual void setExactSolution(AbstractFunctor* exactSolution) override;

    virtual void appendNormsToFile(const double& time, BlockVector solution,
                                   std::ofstream& outFile) override {toBeImplemented();};

    virtual void appendErrorsToFile(const double& time, BlockVector solution,
                                    std::ofstream& outFile) override {toBeImplemented();};

    virtual void printMeshSize(std::string filename) override {toBeImplemented();};

    virtual BlockVector getInitialCondition() override {toBeImplemented();};

protected:
    void assembleConstantMatrices();

    virtual void assembleStiffnessMatrix();

    void assembleDivergenceMatrix();

    virtual void assembleMassMatrixPressure();

    void assembleConvectiveMatrix();

    void assembleJacobianConvectiveMatrix();

    void assembleForcingTerm();

    void assembleForcingTermTimeDerivative();

    void applyNeumannBCs(VectorPtr vector, std::function<double(double)> law);

    void applyNeumannBCsWithExactFunction(VectorPtr vector,
                                          FunctionNeumannType* exactFunction);

    BoundaryConditionPtr createBCHandler(std::function<double(double)> law);

    BoundaryConditionPtr createBCHandlerWithExactFunction(FunctionType*
                                                          exactFunction);

    void computeAndExportErrors(const double& time, BlockVector solutions);

    static double fZero(const double& t,
                        const double& x,
                        const double& y,
                        const double& z,
                        const unsigned int& i);

    static double fOne(const double& t,
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

    static double neumannInflow(const double& t, const double& x, const double& y,
                                const double& z, const unsigned int& i,
                                std::function<double(double)> maxLaw);

    void updateBCs(BoundaryConditionPtr bcToUpdate, FESpacePtr fespace);

    FESpacePtr                              M_velocityFESpace;
    FESpacePtr                              M_pressureFESpace;
    ETFESpaceVelocityPtr                    M_velocityFESpaceETA;
    ETFESpacePressurePtr                    M_pressureFESpaceETA;

    MatrixPtr                               M_A;
    MatrixPtr                               M_B;
    MatrixPtr                               M_Bt;
    MatrixPtr                               M_M;
    MatrixPtr                               M_Mp;
    MatrixPtr                               M_C;
    MatrixPtr                               M_J;
    BlockVector                             M_prevSolution;
    double                                  M_time;

    std::function<double(double)>           M_inflowLaw;
    std::function<double(double)>           M_inflowLawDt;

    ExporterPtr                             M_exporter;
    VectorPtr                               M_velocityExporter;
    VectorPtr                               M_pressureExporter;
    VectorPtr                               M_velocityErrorExporter;
    VectorPtr                               M_pressureErrorExporter;
    VectorPtr                               M_lagrangeMultiplierExporter;

    bool                                    M_useStabilization;
    std::shared_ptr<SUPGStabilization>      M_stabilization;
    bool                                    M_addNoslipBC;
};

}  // namespace RedMA

#endif  // NAVIERSTOKESASSEMBLER_HPP
