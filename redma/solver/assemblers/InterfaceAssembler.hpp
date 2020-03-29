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

#include <redma/solver/assemblers/aAssembler.hpp>
#include <redma/solver/array/BlockMatrix.hpp>
#include <redma/solver/array/BlockVector.hpp>
#include <redma/solver/array/MatrixEp.hpp>
#include <redma/solver/array/VectorEp.hpp>
#include <redma/utils/Exception.hpp>
#include <redma/solver/coupling/BasisFunctionFactory.hpp>

#include <redma/reduced_basis/RBBases.hpp>

#include <lifev/eta/expression/Integrate.hpp>

#define THRESHOLDSTAB       1e-15

namespace RedMA
{

template <class InVectorType, class InMatrixType>
class Interface
{
    typedef aAssembler<InVectorType COMMA InMatrixType>         AssemblerType;
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

template <class InVectorType, class InMatrixType>
class InterfaceAssembler
{
    typedef aAssembler<InVectorType COMMA InMatrixType>         AssemblerType;

public:
    InterfaceAssembler(const DataContainer& data);

    InterfaceAssembler(const DataContainer& data,
                       const Interface<InVectorType, InMatrixType>& interface);

    void setup();

    void buildCouplingMatrices();

    void buildCouplingMatrices(SHP(AssemblerType) assembler,
                               const GeometricFace& face,
                               BlockMatrix<MatrixEp>& matrixT,
                               BlockMatrix<MatrixEp>& matrix);

    void buildStabilizationMatrix(SHP(AssemblerType) assembler,
                                  const GeometricFace& face,
                                  BlockMatrix<InMatrixType>& matrix);

    void buildMapLagrange(SHP(BasisFunctionFunctor) bfs);

    virtual void addContributionRhs(const double& time,
                                    BlockVector<BlockVector<InVectorType>>& rhs,
                                    const BlockVector<BlockVector<InVectorType>>& sol,
                                    const unsigned int& nPrimalBlocks);

    virtual void addContributionJacobianRhs(const double& time,
                                            BlockMatrix<BlockMatrix<InMatrixType>>& jac,
                                            const BlockVector<BlockVector<InVectorType>>& sol,
                                            const unsigned int& nPrimalBlocks);

    inline Interface<InVectorType, InMatrixType> getInterface() const {return M_interface;};

    BlockVector<InVectorType> getZeroVector() const;

    double checkStabilizationTerm(const BlockVector<BlockVector<InVectorType>>& sol,
                                  const unsigned int& nPrimalBlocks);

    inline BlockMatrix<InMatrixType> getFatherBT() const {return M_fatherBT;}

    inline BlockMatrix<InMatrixType> getFatherB() const {return M_fatherB;}

    inline BlockMatrix<InMatrixType> getChildBT() const {return M_childBT;}

    inline BlockMatrix<InMatrixType> getChildB() const {return M_childB;}

protected:
    SHP(LifeV::QuadratureRule) generateQuadratureRule(std::string tag) const;

    std::vector<VectorEp> buildCouplingVectors(SHP(BasisFunctionFunctor) bfs,
                                               const GeometricFace& face,
                                               SHP(aAssembler<InVectorType COMMA InMatrixType>) assembler) const;

    std::vector<InVectorType> buildStabilizationVectorsVelocity(SHP(BasisFunctionFunctor) bfs,
                                                            const GeometricFace& face,
                                                            SHP(aAssembler<InVectorType COMMA InMatrixType>) assembler) const;

    std::vector<InVectorType> buildStabilizationVectorsPressure(SHP(BasisFunctionFunctor) bfs,
                                                            const GeometricFace& face,
                                                            SHP(aAssembler<InVectorType COMMA InMatrixType>) assembler) const;

    std::vector<InVectorType> buildStabilizationVectorsLagrange() const;

    Interface<InVectorType, InMatrixType>           M_interface;

    BlockMatrix<InMatrixType>                       M_identity;
    BlockMatrix<InMatrixType>                       M_fatherBT;
    BlockMatrix<InMatrixType>                       M_fatherB;
    BlockMatrix<InMatrixType>                       M_childBT;
    BlockMatrix<InMatrixType>                       M_childB;
    // this is required in the RB setting to impose weakly dirichlet conditions
    BlockMatrix<MatrixEp>                           M_childBEp;
    BlockMatrix<InMatrixType>                       M_stabChild;
    BlockMatrix<InMatrixType>                       M_stabFather;
    DataContainer                                   M_data;
    SHP(LifeV::MapEpetra)                           M_mapLagrange;
    double                                          M_stabilizationCoupling;
    bool                                            M_isInlet;
};

}

#include "InterfaceAssembler_imp.hpp"

#endif // INTERFACEASSEMBLER_HPP
