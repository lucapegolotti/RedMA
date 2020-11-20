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

    Interface(shp<AssemblerType> assemblerFather, const int& indexFather,
              shp<AssemblerType> assemblerChild, const int& indexChild,
              const unsigned int& interfaceID);

    shp<AssemblerType>      M_assemblerFather;
    shp<AssemblerType>      M_assemblerChild;
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

    void buildCouplingMatrices(shp<AssemblerType> assembler,
                               const GeometricFace& face,
                               shp<BlockMatrix> matrixT,
                               shp<BlockMatrix> matrix);

    void buildStabilizationMatrix(shp<AssemblerType> assembler,
                                  const GeometricFace& face,
                                  shp<BlockMatrix> matrix);

    void buildMapLagrange(shp<BasisFunctionFunctor> bfs);

    virtual void addContributionRhs(const double& time,
                                    shp<BlockVector> rhs,
                                    shp<BlockVector> sol,
                                    const unsigned int& nPrimalBlocks);

    virtual void addContributionJacobianRhs(const double& time,
                                            shp<BlockMatrix> jac,
                                            shp<BlockVector> sol,
                                            const unsigned int& nPrimalBlocks);

    inline Interface getInterface() const {return M_interface;};

    shp<BlockVector> getZeroVector() const;

    double checkStabilizationTerm(const shp<BlockVector>& sol,
                                  const unsigned int& nPrimalBlocks);

    inline shp<BlockMatrix> getFatherBT() const {return M_fatherBT;}

    inline shp<BlockMatrix> getFatherB() const {return M_fatherB;}

    inline shp<BlockMatrix> getChildBT() const {return M_childBT;}

    inline shp<BlockMatrix> getChildB() const {return M_childB;}

protected:
    shp<LifeV::QuadratureRule> generateQuadratureRule(std::string tag) const;

    std::vector<shp<DistributedVector>> buildCouplingVectors(shp<BasisFunctionFunctor> bfs,
                                                             const GeometricFace& face,
                                                             shp<aAssembler> assembler) const;

    std::vector<shp<DistributedVector>> buildStabilizationVectorsVelocity(shp<BasisFunctionFunctor> bfs,
                                                                          const GeometricFace& face,
                                                                          shp<aAssembler> assembler) const;

    std::vector<shp<DistributedVector>> buildStabilizationVectorsPressure(shp<BasisFunctionFunctor> bfs,
                                                                          const GeometricFace& face,
                                                                          shp<aAssembler> assembler) const;

    std::vector<shp<DistributedVector>> buildStabilizationVectorsLagrange() const;

    Interface                              M_interface;
    shp<BlockMatrix>                       M_identity;
    shp<BlockMatrix>                       M_fatherBT;
    shp<BlockMatrix>                       M_fatherB;
    shp<BlockMatrix>                       M_childBT;
    shp<BlockMatrix>                       M_childB;
    // shp<BlockMatrix)                       M_fatherBTreduced;
    // shp<BlockMatrix)                       M_fatherBreduced;
    shp<BlockMatrix>                       M_childBTfe;
    shp<BlockMatrix>                       M_childBfe;
    // this is required in the RB setting to impose weakly dirichlet conditions
    shp<BlockMatrix>                       M_childBEp;
    shp<BlockMatrix>                       M_stabChild;
    shp<BlockMatrix>                       M_stabFather;
    DataContainer                          M_data;
    shp<const LifeV::MapEpetra>            M_mapLagrange;
    double                                 M_stabilizationCoupling;
    bool                                   M_isInlet;
};

}

#endif // INTERFACEASSEMBLER_HPP
