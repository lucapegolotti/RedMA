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

#ifndef GLOBALASSEMBLER_HPP
#define GLOBALASSEMBLER_HPP

#include <fstream>

#include <AbstractFunctor.hpp>

#include <AssemblersFactory.hpp>

#include <TreeStructure.hpp>
#include <AbstractAssembler.hpp>
#include <BlockMatrix.hpp>

#include <lifev/core/array/MapVector.hpp>
#include <lifev/core/array/MatrixEpetraStructured.hpp>
#include <lifev/core/array/VectorEpetraStructured.hpp>
#include <lifev/core/array/MatrixEpetraStructuredUtility.hpp>

namespace RedMA
{

class GlobalAssembler : public AbstractAssembler
{
    typedef std::shared_ptr<TreeNode>                       TreeNodePtr;
    typedef LifeV::MapEpetra                                MapEpetra;
	typedef std::shared_ptr<MapEpetra>                      MapEpetraPtr;
    typedef LifeV::VectorEpetra                             Vector;
    typedef std::shared_ptr<Vector>                         VectorPtr;
    typedef LifeV::MatrixEpetra<double>                     Matrix;
    typedef std::shared_ptr<Matrix>                         MatrixPtr;
    typedef std::shared_ptr<Epetra_Comm>                    commPtr_Type;
    typedef std::shared_ptr<AbstractAssembler>              AbstractAssemblerPtr;

    typedef std::function<double(double const&,
                                 double const&,
                                 double const&,
                                 double const&,
                                 unsigned int const& )>     Function;

    typedef std::shared_ptr<AbstractMatrix>                AbstractMatrixPtr;
    typedef std::shared_ptr<AbstractVector>                AbstractVectorPtr;
    typedef std::shared_ptr<BlockVector>                   BlockVectorPtr;

public:
    GlobalAssembler(const GetPot& datafile, commPtr_Type comm,
                    bool verbose = false);

    virtual void setup() override;

    virtual unsigned int numberOfBlocks() override {return M_dimensionsVector.size();};

    // this does not make sense for global assembler
    virtual unsigned int numberOfComponents() override;

    virtual void applyBCsMatrix(MatrixPtr matrix, const double& diagonalCoefficient,
                                const unsigned int& iblock, const unsigned int& jblock) override {};

    virtual MatrixPtr getUpdateMass(const unsigned int& blockrow,
                                    const unsigned int& blockcol) override {};

    virtual MatrixPtr getUpdateMassJac(const unsigned int& blockrow,
                                       const unsigned int& blockcol) override {};

    virtual MatrixPtr getUpdateMassJacVelocity(const unsigned int& blockrow,
                                               const unsigned int& blockcol) override {};

    virtual std::vector<VectorPtr> initialCondition() override {};

    virtual MatrixPtr getJacobian(const unsigned int& blockrow,
                                  const unsigned int& blockcol) override {};

    virtual MatrixPtr getJacobianPrec(const unsigned int& blockrow,
                                    const unsigned int& blockcol) override {};

    virtual void setTimeAndPrevSolution(const double& time, BlockVector solution,
                                        bool doAssembly = true) override;

    virtual void applyBCsRhsRosenbrock(std::vector<VectorPtr> rhs,
                                       std::vector<VectorPtr> utilde,
                                       const double& time,
                                       const double& dt,
                                       const double& alphai,
                                       const double& gammai) override {};

    virtual std::vector<VectorPtr> computeF() override {};

    virtual std::vector<VectorPtr> computeFder() override {};

    virtual std::vector<double> computeNorms(std::vector<VectorPtr> solutions) override {};

    virtual std::vector<double> computeErrors(std::vector<VectorPtr> solutions,
                                              const double& time) override {};

    virtual void exportSolutions(const double& time, std::vector<VectorPtr> solutions) override {};

    virtual void setTreeStructure(TreeStructure& tree) override;

    void buildPrimalStructures(TreeStructure& tree);

    void buildDualStructures(TreeStructure& tree);

    MapEpetraPtr getGlobalMap() const override;

    virtual AbstractMatrixPtr getMassMatrix() override;

    BlockMatrix getGlobalMassJac();

    BlockMatrix getGlobalMassJacVelocity();

    BlockMatrix getJacobianF(bool addCoupling,
                                   double* diagonalCoefficient = nullptr);

    BlockMatrix getJacobianFprec(bool addCoupling,
                                       double* diagonalCoefficient = nullptr);

    VectorPtr computeF_();

    VectorPtr computeFder_();

    // the diagonal coefficient is for the boundary conditions (if null, no
    // bcs are applied)
    AbstractMatrixPtr assembleMassMatrix(double* diagonalCoefficient = nullptr) override;

    void applyBCsRhsRosenbrock(VectorPtr rhs, VectorPtr utilde,
                               const double& time, const double& dt,
                               const double& alphai, const double& gammai);

    virtual void applyBCsBackwardEuler(BlockVector rhs, const double& coeff,
                                       const double& time) override;

    virtual void setLawInflow(std::function<double(double)> maxLaw) override;

    virtual void setLawDtInflow(std::function<double(double)> maxLawDt) override;

    virtual void setExactSolution(AbstractFunctor* exactSolution) override;

    void exportSolutions(const double& time, VectorPtr solution);

    virtual void appendNormsToFile(const double& time, BlockVector solution,
                                   std::ofstream& outFile) override;

    virtual void appendErrorsToFile(const double& time, BlockVector solution,
                                    std::ofstream& outFile) override;

    virtual void setTimeIntegrationOrder(unsigned int order) override;

    virtual void setTimestep(double dt) override;

    virtual void postProcess() override;

    BlockVector getInitialCondition() override;

    void printMeshSize(std::string filename) override;

    void setForcingFunction(Function forcingFunction,
                            Function forcingFunctionDt);

private:
    template<typename FunctionType>
    void fillGlobalMatrix(BlockMatrix& matrixToFill,
                          bool addCoupling,
                          FunctionType getMatrixMethod,
                          double* diagonalCoefficient);

    template<typename FunctionType>
    void fillGlobalVector(VectorPtr& vectorToFill,
                          FunctionType getVectorMethod);

    std::vector<std::pair<unsigned int, AbstractAssemblerPtr> > M_assemblersVector;
    std::map<unsigned int, AbstractAssemblerPtr>                M_assemblersMap;
    GetPot                                                      M_datafile;
    commPtr_Type                                                M_comm;
    bool                                                        M_verbose;
    MapEpetraPtr                                                M_globalMap;
    std::shared_ptr<BlockMatrix>                                M_massMatrix;
    std::shared_ptr<BlockMatrix>                                M_updatedMassMatrix;
    std::vector<unsigned int>                                   M_dimensionsVector;
    std::vector<std::pair<unsigned int, unsigned int> >         M_interfaces;
    std::vector<unsigned int>                                   M_offsets;
    unsigned int                                                M_nPrimalBlocks;
    std::vector<MapEpetraPtr>                                   M_maps;
    bool                                                        M_updateMassMatrix;
    TreeStructure                                               M_tree;
};

template<typename FunctionType>
void
GlobalAssembler::
fillGlobalMatrix(BlockMatrix& matrixToFill, bool addCoupling,
                 FunctionType getMatrixMethod, double* diagonalCoefficient)
{
    toBeImplemented();

    // typedef std::vector<std::pair<unsigned int, AbstractAssemblerPtr> >
    //             AssemblersVector;
    //
    // unsigned int totalNumberBlocks = M_nPrimalBlocks + M_interfaces.size();
    // matrixToFill.resize(totalNumberBlocks, totalNumberBlocks);
    // matrixToFill.setMaps(M_maps, M_maps);
    // // we start with the primal blocks
    // unsigned int countBlocks = 0;
    // for (typename AssemblersVector::iterator it = M_assemblersVector.begin();
    //      it != M_assemblersVector.end(); it++)
    // {
    //     unsigned int blockIndex = it->first;
    //     unsigned int numberBlocks = it->second->numberOfBlocks();
    //     for (int i = 0; i < numberBlocks; i++)
    //     {
    //         for (int j = 0; j < numberBlocks; j++)
    //         {
    //             AbstractAssembler& curAssembler = *it->second;
    //             MatrixPtr localMatrix = (curAssembler.*getMatrixMethod)(i, j);
    //             if (diagonalCoefficient)
    //                 curAssembler.applyBCsMatrix(localMatrix, *diagonalCoefficient,
    //                                             i, j);
    //             // Attention: this does not work if number of blocks is not constant
    //             // over all the domains
    //             matrixToFill.copyBlock(countBlocks + i,
    //                                    countBlocks + j,
    //                                    localMatrix);
    //         }
    //     }
    //     countBlocks += numberBlocks;
    // }
    //
    // // if required add the coupling blocks
    // if (addCoupling)
    // {
    //     typedef std::vector<std::pair<unsigned int, unsigned int> >
    //             InterfacesVector;
    //
    //     unsigned int offset = M_nPrimalBlocks;
    //
    //     for (InterfacesVector::iterator it = M_interfaces.begin();
    //          it != M_interfaces.end(); it++)
    //     {
    //         unsigned int indices[2] = {it->first, it->second};
    //
    //         for (int i = 0; i < 2; i++)
    //         {
    //             AbstractAssembler& curAssembler = *M_assemblersMap[indices[i]];
    //             MatrixPtr Qt = curAssembler.getQT(indices[(i+1) % 2]);
    //             unsigned int blockCoupling = curAssembler.getIndexCoupling();
    //
    //             unsigned int blockIndex = blockCoupling +
    //                             indices[i] * curAssembler.numberOfBlocks();
    //
    //             curAssembler.applyBCsMatrix(Qt, 0.0,
    //                                         blockCoupling, blockCoupling);
    //             matrixToFill.copyBlock(blockIndex, offset, Qt);
    //
    //             MatrixPtr Q = curAssembler.getQ(indices[(i+1) % 2]);
    //             matrixToFill.copyBlock(offset, blockIndex, Q);
    //         }
    //         offset++;
    //     }
    // }
}

template<typename FunctionType>
void
GlobalAssembler::
fillGlobalVector(VectorPtr& vectorToFill, FunctionType getVectorMethod)
{
    using namespace LifeV::MatrixEpetraStructuredUtility;

    typedef std::vector<std::pair<unsigned int, AbstractAssemblerPtr> >
                AssemblersVector;

    typedef std::vector<MapEpetraPtr>                    MapVector;

    vectorToFill.reset(new Vector(*M_globalMap));
    vectorToFill->zero();

    unsigned int offset = 0;
    for (typename AssemblersVector::iterator it = M_assemblersVector.begin();
         it != M_assemblersVector.end(); it++)
    {
        std::vector<VectorPtr> localSolutions;
        MapVector maps = it->second->getPrimalMapVector();

        AbstractAssembler& curAssembler = *it->second;
        std::vector<VectorPtr> localVectors;
        localVectors = (curAssembler.*getVectorMethod)();

        unsigned int index = 0;

        for (MapVector::iterator itmap = maps.begin();
             itmap != maps.end(); itmap++)
        {
            LifeV::MapEpetra& curLocalMap = **itmap;
            // we fill only if the vector corresponding to map i exists!
            if (localVectors[index])
            {
                vectorToFill->subset(*localVectors[index], curLocalMap, 0,
                                     offset);
            }
            offset += curLocalMap.mapSize();
            index++;
        }

        // deal with the dual part. This is trickier because more assemblers
        // contribute to the same subsets of vectorfill
        index = 0;
        maps = it->second->getDualMapVector();
        std::vector<unsigned int> indices = it->second->getInterfacesIndices();
        for (MapVector::iterator itmap = maps.begin();
             itmap != maps.end(); itmap++)
        {
            LifeV::MapEpetra& curLocalMap = **itmap;
            VectorPtr aux(new Vector(curLocalMap));
            aux->subset(*vectorToFill, curLocalMap,
                        M_offsets[M_nPrimalBlocks + indices[index]], 0);
            *aux += 0;
            localVectors[index + it->second->numberOfBlocks()]->zero();
            vectorToFill->subset(*aux, curLocalMap, 0,
                                  M_offsets[M_nPrimalBlocks + indices[index]]);
            index++;
        }
    }
}

}  // namespace RedMA

#endif  // GLOBALASSEMBLER_HPP
