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

#ifndef INTERFACEASSEMBLER_HPP
#define INTERFACEASSEMBLER_HPP

#include <redma/RedMA.hpp>

#include <redma/assemblers/abstract/aAssembler.hpp>
#include <redma/array/BlockMatrix.hpp>
#include <redma/array/BlockVector.hpp>
#include <redma/array/SparseMatrix.hpp>
#include <redma/array/DistributedVector.hpp>
#include <redma/utils/Exception.hpp>
#include <redma/coupling_basis_functions/BasisFunctionFactory.hpp>

#include <redma/reduced_basis/RBBases.hpp>

#include <lifev/eta/expression/Integrate.hpp>

#define THRESHOLDSTAB       1e-15

namespace RedMA
{

class Interface
{
    typedef aAssembler         AssemblerType;
public:
    Interface ();

    Interface(SHP(AssemblerType) assemblerFather, const int& indexFather,
              SHP(AssemblerType) assemblerChild, const int& indexChild,
              const unsigned int& interfaceID);

    SHP(AssemblerType)      M_assemblerFather;
    SHP(AssemblerType)      M_assemblerChild;
    int                     M_indexFather;
    int                     M_indexChild;
    unsigned int            M_ID;
    unsigned int            M_indexOutlet;
};

class InterfaceAssembler
{
    typedef aAssembler         AssemblerType;

public:

    InterfaceAssembler(const DataContainer& data);

    InterfaceAssembler(const DataContainer& data,
                       const Interface& interface);

    virtual ~InterfaceAssembler() {}

    void setup();

    void buildCouplingMatrices();

    void buildCouplingMatrices(SHP(AssemblerType) assembler,
                               const GeometricFace& face,
                               SHP(BlockMatrix) matrixT,
                               SHP(BlockMatrix) matrix);

    void buildStabilizationMatrix(SHP(AssemblerType) assembler,
                                  const GeometricFace& face,
                                  SHP(BlockMatrix) matrix);

    void buildMapLagrange(SHP(BasisFunctionFunctor) bfs);

    virtual void addContributionRhs(const double& time,
                                    SHP(BlockVector) rhs,
                                    SHP(BlockVector) sol,
                                    const unsigned int& nPrimalBlocks);

    virtual void addContributionJacobianRhs(const double& time,
                                            SHP(BlockMatrix) jac,
                                            SHP(BlockVector) sol,
                                            const unsigned int& nPrimalBlocks);

    inline Interface getInterface() const {return M_interface;};

    SHP(BlockVector) getZeroVector() const;

    double checkStabilizationTerm(const SHP(BlockVector)& sol,
                                  const unsigned int& nPrimalBlocks);

    inline SHP(BlockMatrix) getFatherBT() const {return M_fatherBT;}

    inline SHP(BlockMatrix) getFatherB() const {return M_fatherB;}

    inline SHP(BlockMatrix) getChildBT() const {return M_childBT;}

    inline SHP(BlockMatrix) getChildB() const {return M_childB;}

protected:
    SHP(LifeV::QuadratureRule) generateQuadratureRule(std::string tag) const;

    std::vector<SHP(DistributedVector)> buildCouplingVectors(SHP(BasisFunctionFunctor) bfs,
                                                             const GeometricFace& face,
                                                             SHP(aAssembler) assembler) const;

    std::vector<SHP(DistributedVector)> buildStabilizationVectorsVelocity(SHP(BasisFunctionFunctor) bfs,
                                                                          const GeometricFace& face,
                                                                          SHP(aAssembler) assembler) const;

    std::vector<SHP(DistributedVector)> buildStabilizationVectorsPressure(SHP(BasisFunctionFunctor) bfs,
                                                                          const GeometricFace& face,
                                                                          SHP(aAssembler) assembler) const;

    std::vector<SHP(DistributedVector)> buildStabilizationVectorsLagrange() const;

    Interface                              M_interface;
    SHP(BlockMatrix)                       M_identity;
    SHP(BlockMatrix)                       M_fatherBT;
    SHP(BlockMatrix)                       M_fatherB;
    SHP(BlockMatrix)                       M_childBT;
    SHP(BlockMatrix)                       M_childB;
    // this is required in the RB setting to impose weakly dirichlet conditions
    SHP(BlockMatrix)                       M_childBEp;
    SHP(BlockMatrix)                       M_stabChild;
    SHP(BlockMatrix)                       M_stabFather;
    DataContainer                          M_data;
    SHP(LifeV::MapEpetra)                  M_mapLagrange;
    double                                 M_stabilizationCoupling;
    bool                                   M_isInlet;
};

}

#endif // INTERFACEASSEMBLER_HPP
