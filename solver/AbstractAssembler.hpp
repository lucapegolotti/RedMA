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
#include <FourierBasisFunction.hpp>

#include <lifev/core/array/MapEpetra.hpp>
#include <lifev/core/filter/GetPot.hpp>
#include <lifev/core/array/MatrixEpetra.hpp>

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
    typedef LifeV::MatrixSmall<3, 3>                        Matrix3D;
    typedef LifeV::FESpace<Mesh, MapEpetra>                FESpace;
    typedef std::shared_ptr<FESpace>                       FESpacePtr;
    typedef LifeV::VectorEpetra                            Vector;
    typedef std::shared_ptr<Vector>                        VectorPtr;
    typedef LifeV::MatrixEpetra<double>                    Matrix;
    typedef std::shared_ptr<Matrix>                        MatrixPtr;
    typedef LifeV::ETFESpace<Mesh, MapEpetra, 3, 1>        ETFESpaceCoupling;
    typedef std::shared_ptr<ETFESpaceCoupling>             ETFESpaceCouplingPtr;

public:
    AbstractAssembler(const GetPot& datafile, commPtr_Type comm,
                      const TreeNodePtr& treeNode, bool verbose = false);

    void addPrimalMaps(MapEpetraPtr& globalMap,
                       std::vector<unsigned int>& dimensions);

    // this method should be used to:
    // 1) create the finite element spaces
    // 2) assemble the constant matrices
    virtual void setup() = 0;

    virtual unsigned int numberOfBlocks() = 0;

    // number of components of the variable involved in the coupling
    virtual unsigned int numberOfComponents() = 0;

    std::vector<MapEpetraPtr> getPrimalMapVector();

    std::vector<MapEpetraPtr> getDualMapVector();

    void buildLagrangeMultiplierBasisFourier(const unsigned int& frequencies,
                                             const unsigned int& nComponents,
                                             ETFESpaceCouplingPtr couplingFespace,
                                             MapEpetraPtr primalMap,
                                             FESpacePtr primalFespace,
                                             GeometricFace face,
                                             const unsigned int& faceFlag);

    // this function must be called on the father
    void assembleCouplingMatrices(AbstractAssembler& child,
                                  const unsigned int& indexOutlet,
                                  const unsigned int& interfaceIndex,
                                  MapEpetraPtr& globalMap,
                                  std::vector<unsigned int>& dimensions);

    MatrixPtr getQT(const unsigned int& flag);

    MatrixPtr getQ(const unsigned int& flag);

    unsigned int getIndexCoupling();

    std::vector<unsigned int> getInterfacesIndices();

private:
    static void gramSchmidt(VectorPtr* basis1, MatrixPtr massMatrix1,
                            VectorPtr* basis2, MatrixPtr massMatrix2,
                            unsigned int& nVectors);

    static double dotProd(VectorPtr* basis1, VectorPtr* basis2, unsigned int index1,
                          unsigned int index2, MatrixPtr mass1, MatrixPtr mass2);

    VectorPtr* assembleCouplingVectorsFourier(const unsigned int& frequencies,
                                              const unsigned int& nBasisFunctions,
                                              GeometricFace face,
                                              const double& coeff);

    void fillMatricesWithVectors(VectorPtr* couplingVectors,
                                 const unsigned int& nBasisFunctions,
                                 MapEpetraPtr lagrangeMap,
                                 const unsigned int& flagAdjacentDomain);

    MatrixPtr assembleBoundaryMatrix(GeometricFace face);

protected:
    TreeNodePtr                         M_treeNode;
    std::vector<MapEpetraPtr>           M_primalMaps;
    std::vector<MapEpetraPtr>           M_dualMaps;
    GetPot                              M_datafile;
    commPtr_Type                        M_comm;
    bool                                M_verbose;
    // maps of the coupling matrices (key = flag of corresponding face)
    std::map<unsigned int, MatrixPtr>   M_mapQTs;
    std::map<unsigned int, MatrixPtr>   M_mapQs;
    ETFESpaceCouplingPtr                M_couplingFESpaceETA;
    // index of the block to which the coupling must be applied
    unsigned int                        M_indexCoupling;
    std::vector<unsigned int>           M_interfacesIndices;
};

}  // namespace RedMA

#endif  // ABSTRACTASSEMBLER_HPP
