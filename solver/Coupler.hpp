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

#ifndef COUPLER_HPP
#define COUPLER_HPP

#include <TreeStructure.hpp>
#include <PrintLog.hpp>
#include <Exception.hpp>
#include <FourierBasisFunction.hpp>
#include <ZernikeBasisFunction.hpp>

#include <lifev/core/array/MapEpetra.hpp>
#include <lifev/core/filter/GetPot.hpp>
#include <lifev/core/array/MatrixEpetra.hpp>
#include <lifev/core/interpolation/Interpolation.hpp>
#include <lifev/core/mesh/MeshUtility.hpp>

#include <lifev/core/fem/BCHandler.hpp>

#include <Epetra_LAPACK.h>

namespace RedMA
{

class Coupler
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
    typedef LifeV::ETFESpace<Mesh, MapEpetra, 3, 1>        ETFESpaceCoupling;
    typedef std::shared_ptr<ETFESpaceCoupling>             ETFESpaceCouplingPtr;
    typedef std::shared_ptr<LifeV::BCHandler>              BoundaryConditionPtr;
    typedef std::function<double(double const&,
                                 double const&,
                                 double const&,
                                 double const&,
                                 unsigned int const& )>    Function;
    typedef LifeV::Interpolation                           Interpolation;
    typedef std::shared_ptr<Interpolation>                 InterpolationPtr;
public:
    Coupler(commPtr_Type comm, bool verbose);

    MatrixPtr getQT(const unsigned int& flag);

    MatrixPtr getQ(const unsigned int& flag);

    static void gramSchmidt(VectorPtr* basis1, VectorPtr* basis2,
                            MatrixPtr massMatrix1, MatrixPtr massMatrix2,
                            unsigned int& nVectors);

    void POD(VectorPtr*& basis1, VectorPtr*& basis2,
             MatrixPtr massMatrix1, MatrixPtr massMatrix2,
             unsigned int& nVectors,
             double tol,
             unsigned int offset = 0);

    static double dotProd(VectorPtr* basis1, VectorPtr* basis2, unsigned int index1,
                          unsigned int index2, MatrixPtr mass1 = nullptr,
                          MatrixPtr mass2 = nullptr);

    static double* computeCorrelationMatrix(VectorPtr* basis1, VectorPtr* basis2,
                                            MatrixPtr mass1, MatrixPtr mass2,
                                            const unsigned int& nBasisFunctions);

    VectorPtr* assembleCouplingVectors(std::shared_ptr<BasisFunctionFunctor> bf,
                                       GeometricFace face,
                                       const double& coeff,
                                       FESpacePtr couplingFespace,
                                       ETFESpaceCouplingPtr couplingFESpaceETA,
                                       VectorPtr* otherInterfaceVectors = nullptr,
                                       InterpolationPtr interpolator = nullptr);

    void fillMatricesWithVectors(VectorPtr* couplingVectors,
                                 const unsigned int& nBasisFunctions,
                                 MapEpetraPtr lagrangeMap,
                                 MapEpetraPtr map,
                                 unsigned int numberOfComponents,
                                 const unsigned int& flagAdjacentDomain);

    void fillMatrixWithVectorsInterpolated(VectorPtr* couplingVectors,
                                           const unsigned int& nBasisFunctions,
                                           MapEpetraPtr lagrangeMap,
                                           MapEpetraPtr map,
                                           unsigned int numberOfComponents,
                                           const unsigned int& flagAdjacentDomain);

    std::map<unsigned int, MatrixPtr>& getMapsQTs();

    std::map<unsigned int, MatrixPtr>& getMapsQTsInterpolated();

    std::map<unsigned int, MatrixPtr>& getMapsQs();

protected:
    commPtr_Type                        M_comm;
    bool                                M_verbose;
    // maps of the coupling matrices (key = flag of corresponding face)
    std::map<unsigned int, MatrixPtr>   M_mapQTs;
    std::map<unsigned int, MatrixPtr>   M_mapQs;
    std::map<unsigned int, MatrixPtr>   M_mapQTsInterpolated;
};

}  // namespace RedMA

#endif  // COUPLER_HPP
