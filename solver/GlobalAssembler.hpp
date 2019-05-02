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

#include <TreeStructure.hpp>

#include <lifev/core/array/MapVector.hpp>
#include <lifev/core/array/MatrixEpetraStructured.hpp>
#include <lifev/core/array/VectorEpetraStructured.hpp>
#include <lifev/core/array/MatrixEpetraStructuredUtility.hpp>

namespace RedMA
{

template <class AssemblerType>
class GlobalAssembler
{
    typedef std::shared_ptr<TreeNode>                       TreeNodePtr;
    typedef std::shared_ptr<AssemblerType>                  AssemblerTypePtr;
    typedef LifeV::MapEpetra                                MapEpetra;
	typedef std::shared_ptr<MapEpetra>                      MapEpetraPtr;
    typedef LifeV::VectorEpetra                             Vector;
    typedef std::shared_ptr<Vector>                         VectorPtr;
    typedef LifeV::MatrixEpetra<double>                     Matrix;
    typedef std::shared_ptr<Matrix>                         MatrixPtr;
    typedef std::shared_ptr<Epetra_Comm>                    commPtr_Type;

public:
    GlobalAssembler(const GetPot& datafile, commPtr_Type comm,
                    bool verbose = false);

    void buildPrimalStructures(TreeStructure& tree);

    MapEpetraPtr getGlobalMap() const;

    MatrixPtr getGlobalMass() const;

    MatrixPtr getJacobianF(double* diagonalCoefficient = nullptr);

    VectorPtr computeF();

    VectorPtr computeFder();

    // the diagonal coefficient is for the boundary conditions (if null, no
    // bcs are applied)
    MatrixPtr assembleGlobalMass(double* diagonalCoefficient = nullptr);

    void setTimeAndPrevSolution(const double& time, VectorPtr solution);

    void applyBCsRhsRosenbrock(VectorPtr rhs, VectorPtr utilde,
                               const double& time, const double& dt,
                               const double& alphai, const double& gammai);

    void setMaxVelocityLawInflow(std::function<double(double)> maxLaw);

    void setMaxVelocityDtLawInflow(std::function<double(double)> maxLawDt);

    void exportSolutions(const double& t, VectorPtr solution);

private:
    template<typename FunctionType>
    void fillGlobalMatrix(MatrixPtr& matrixToFill,
                          FunctionType getMatrixMethod,
                          double* diagonalCoefficient);

    template<typename FunctionType>
    void fillGlobalVector(VectorPtr& vectorToFill,
                          FunctionType getVectorMethod);

    std::vector<std::pair<unsigned int, AssemblerTypePtr> > M_assemblersVector;
    GetPot                                                  M_datafile;
    commPtr_Type                                            M_comm;
    bool                                                    M_verbose;
    MapEpetraPtr                                            M_globalMap;
    MatrixPtr                                               M_massMatrix;
    std::vector<unsigned int>                               M_dimensionsVector;
};

}  // namespace RedMA

#include "GlobalAssembler_imp.hpp"

#endif  // GLOBALASSEMBLER_HPP
