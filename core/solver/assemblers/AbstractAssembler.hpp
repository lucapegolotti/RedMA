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
#include <Exception.hpp>
#include <BasisFunctionFactory.hpp>
#include <Coupler.hpp>
#include <AbstractFunctor.hpp>
#include <TreeStructure.hpp>

#include <lifev/core/array/MapEpetra.hpp>
#include <lifev/core/filter/GetPot.hpp>
#include <lifev/core/array/MatrixEpetra.hpp>
#include <lifev/core/interpolation/Interpolation.hpp>
#include <lifev/core/mesh/MeshUtility.hpp>
#include <AbstractMatrix.hpp>
#include <BlockVector.hpp>

#include <lifev/core/fem/BCHandler.hpp>

#include <Epetra_LAPACK.h>

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
    typedef LifeV::MatrixSmall<3, 3>                       Matrix3D;
    typedef LifeV::FESpace<Mesh, MapEpetra>                FESpace;
    typedef std::shared_ptr<FESpace>                       FESpacePtr;
    typedef LifeV::VectorEpetra                            Vector;
    typedef std::shared_ptr<Vector>                        VectorPtr;
    typedef LifeV::MatrixEpetra<double>                    Matrix;
    typedef std::shared_ptr<Matrix>                        MatrixPtr;
    typedef LifeV::ETFESpace<Mesh, MapEpetra, 3, 1>        ETFESpaceCouplingScalar;
    typedef std::shared_ptr<ETFESpaceCouplingScalar>       ETFESpaceCouplingScalarPtr;
    typedef LifeV::ETFESpace<Mesh, MapEpetra, 3, 3>        ETFESpaceCoupling;
    typedef std::shared_ptr<ETFESpaceCoupling>             ETFESpaceCouplingPtr;
    typedef std::shared_ptr<LifeV::BCHandler>              BoundaryConditionPtr;
    typedef std::function<double(double const&,
                                 double const&,
                                 double const&,
                                 double const&,
                                 unsigned int const& )>    Function;
    typedef LifeV::Interpolation                           Interpolation;
    typedef std::shared_ptr<Interpolation>                 InterpolationPtr;
    typedef std::shared_ptr<AbstractMatrix>                AbstractMatrixPtr;
    typedef std::shared_ptr<AbstractVector>                AbstractVectorPtr;
    typedef std::shared_ptr<BlockVector>                   BlockVectorPtr;
public:
    AbstractAssembler(const GetPot& datafile, commPtr_Type comm,
                      const TreeNodePtr& treeNode, bool verbose = false);

    void addPrimalMaps(MapEpetraPtr& globalMap,
                       std::vector<MapEpetraPtr>& maps,
                       std::vector<unsigned int>& dimensions);

    // this method should be used to:
    // 1) create the finite element spaces
    // 2) assemble the constant matrices
    virtual void setup() = 0;

    virtual unsigned int numberOfBlocks() = 0;

    // number of components of the variable involved in the coupling
    virtual unsigned int numberOfComponents() = 0;

    virtual MapEpetraPtr getGlobalMap() const = 0;

    virtual AbstractMatrixPtr assembleMassMatrix(double* diagonalCoefficient = nullptr) = 0;

    virtual void applyBCsMatrix(MatrixPtr matrix, const double& diagonalCoefficient,
                        const unsigned int& iblock, const unsigned int& jblock) = 0;

    virtual MatrixPtr getUpdateMass(const unsigned int& blockrow,
                                    const unsigned int& blockcol) = 0;

    virtual MatrixPtr getUpdateMassJac(const unsigned int& blockrow,
                                       const unsigned int& blockcol) = 0;

    virtual MatrixPtr getUpdateMassJacVelocity(const unsigned int& blockrow,
                                               const unsigned int& blockcol) = 0;

    virtual MatrixPtr getJacobian(const unsigned int& blockrow,
                                  const unsigned int& blockcol) = 0;

    virtual MatrixPtr getJacobianPrec(const unsigned int& blockrow,
                                      const unsigned int& blockcol) = 0;

    virtual AbstractMatrixPtr getMassMatrix() = 0;

    virtual std::vector<VectorPtr> initialCondition() = 0;

    virtual std::vector<VectorPtr> computeF() = 0;

    virtual std::vector<VectorPtr> computeFder() = 0;

    virtual void setTimeAndPrevSolution(const double& time,
                                        BlockVector solution,
                                        bool assembleBlocks = true) = 0;

    virtual void applyBCsRhsRosenbrock(std::vector<VectorPtr> rhs,
                                       std::vector<VectorPtr> utilde,
                                       const double& time,
                                       const double& dt,
                                       const double& alphai,
                                       const double& gammai) = 0;

    virtual void applyBCsBackwardEuler(BlockVector rhs, const double& coeff,
                                       const double& time) = 0;

    virtual void setLawInflow(std::function<double(double)> maxLaw) = 0;

    virtual void setLawDtInflow(std::function<double(double)> maxLawDt) = 0;

    virtual std::vector<double> computeNorms(std::vector<VectorPtr> solutions) = 0;

    virtual std::vector<double> computeErrors(std::vector<VectorPtr> solutions,
                                              const double& time) = 0;

    virtual void exportSolutions(const double& time, std::vector<VectorPtr> solutions) = 0;

    static std::string normFileFirstLine();

    static std::string errorsFileFirstLine();

    virtual void appendNormsToFile(const double& time, BlockVector solution,
                                   std::ofstream& outFile) = 0;

    virtual void appendErrorsToFile(const double& time, BlockVector solution,
                                    std::ofstream& outFile) = 0;

    virtual void printMeshSize(std::string filename) = 0;

    virtual BlockVector getInitialCondition() = 0;

    // this is an ugly way to allow for easy polymorphism between assemblers and
    // global assembler
    virtual void setTreeStructure(TreeStructure& tree) {}

    std::vector<MapEpetraPtr> getPrimalMapVector();

    std::vector<MapEpetraPtr> getDualMapVector();

    // this function must be called on the father
    void assembleCouplingMatrices(AbstractAssembler& child,
                                  const unsigned int& indexOutlet,
                                  const unsigned int& interfaceIndex,
                                  MapEpetraPtr& globalMap,
                                  std::vector<MapEpetraPtr>& maps,
                                  std::vector<unsigned int>& dimensions);

    MatrixPtr getQT(const unsigned int& flag);

    MatrixPtr getQ(const unsigned int& flag);

    std::map<unsigned int, MatrixPtr>& getMapsQs();

    std::map<unsigned int, MatrixPtr>& getMapsQTs();

    unsigned int getIndexCoupling();

    std::vector<unsigned int> getInterfacesIndices();

    virtual void setTimeIntegrationOrder(unsigned int order);

    virtual void setTimestep(double dt);

    virtual void postProcess(){};

    VectorPtr reconstructLagrangeMultipliers(std::vector<VectorPtr> solutions, unsigned int offset);

    inline FESpacePtr getCouplingFESpace() {return M_couplingFESpace;};

    virtual void setExactSolution(AbstractFunctor* exactSolution);

    double getMeshSize() {return M_meshSize;};

    void setForcingFunction(Function forcingFunction, Function functionDt);

    void toBeImplemented() const {throw new Exception("This function must still be implemented!");};

private:

    MatrixPtr assembleBoundaryMatrix(GeometricFace face);

    void multiplyVectorsByMassMatrix(VectorPtr* couplingVectors,
                                     const unsigned int& nBasisFunctions,
                                     MatrixPtr massMatrix);

    std::shared_ptr<BasisFunctionFunctor> castBasisFunction(GeometricFace inlet);

    void assembleCouplingMatricesInterpolation(AbstractAssembler& child,
                                               const unsigned int& indexOutlet,
                                               const unsigned int& interfaceIndex,
                                               std::shared_ptr<BasisFunctionFunctor> basisFunction,
                                               MapEpetraPtr& globalMap,
                                               std::vector<MapEpetraPtr>& maps,
                                               std::vector<unsigned int>& dimensions);

    MapEpetraPtr buildLagrangeMultiplierMap(const unsigned int nBasisFunctions,
                                            AbstractAssembler& child,
                                            MapEpetraPtr& globalMap,
                                            std::vector<MapEpetraPtr>& maps,
                                            std::vector<unsigned int>& dimensions);

protected:
    TreeNodePtr                         M_treeNode;
    std::vector<MapEpetraPtr>           M_primalMaps;
    std::vector<MapEpetraPtr>           M_dualMaps;
    GetPot                              M_datafile;
    commPtr_Type                        M_comm;
    bool                                M_verbose;
    FESpacePtr                          M_couplingFESpace;
    ETFESpaceCouplingScalarPtr          M_couplingFESpaceScalarETA;
    ETFESpaceCouplingPtr                M_couplingFESpaceETA;
    // index of the block to which the coupling must be applied
    unsigned int                        M_indexCoupling;
    std::vector<unsigned int>           M_interfacesIndices;
    VectorPtr                           M_couplingVector;
    Coupler                             M_coupler;
    unsigned                            M_timeIntegrationOrder;
    double                              M_dt;
    AbstractFunctor*                    M_exactSolution;
    double                              M_meshSize;
    Function                            M_forceFunction;
    Function                            M_forceTimeDerFunction;
    VectorPtr                           M_forcingTerm;
    VectorPtr                           M_forcingTermTimeDer;
};

}  // namespace RedMA

#endif  // ABSTRACTASSEMBLER_HPP
