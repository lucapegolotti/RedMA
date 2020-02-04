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

#ifndef aASSEMBLER_HPP
#define aASSEMBLER_HPP

#include <redma/RedMA.hpp>
#include <redma/solver/array/BlockMatrix.hpp>
#include <redma/solver/boundary_conditions/BCManager.hpp>
#include <redma/solver/time_marching_algorithms/aFunctionProvider.hpp>
#include <redma/geometry/TreeStructure.hpp>

#include <redma/solver/problem/DataContainer.hpp>

namespace RedMA
{

template <class InVectorType, class InMatrixType>
class aAssembler : public aFunctionProvider<InVectorType, InMatrixType>
{
public:
    aAssembler(const DataContainer& datafile);

    aAssembler(const DataContainer& datafile, SHP(TreeNode) node);

    inline SHP(BCManager) getBCManager() const {return M_bcManager;}

    virtual void setup() = 0;

    virtual void exportSolution(const double& t,
                                const BlockVector<InVectorType>& sol) = 0;

    virtual void postProcess(const BlockVector<InVectorType>& sol) = 0;

    virtual BlockMatrix<InMatrixType> getMass(const double& time,
                                      const BlockVector<InVectorType>& sol) = 0;

    virtual BlockVector<InVectorType> getRightHandSide(const double& time,
                                      const BlockVector<InVectorType>& sol) = 0;

    virtual BlockMatrix<InMatrixType> getJacobianRightHandSide(const double& time,
                                      const BlockVector<InVectorType>& sol) = 0;

    virtual BlockVector<InVectorType> getZeroVector() const = 0;

    // this must be implemented by the inner assemblers
    virtual inline SHP(FESPACE) getFESpaceBCs() const {return nullptr;}

    virtual inline unsigned int getComponentBCs() const {return 0;}

    virtual inline SHP(TreeNode) getTreeNode() const {return M_treeNode;}

    virtual inline unsigned int getNumComponents() const {return M_nComponents;}

    virtual inline EPETRACOMM getComm() const {return M_comm;}

    virtual inline SHP(ETFESPACE3) getETFESpaceCoupling() const {return nullptr;}

    virtual inline SHP(ETFESPACE1) getETFESpaceSecondary() const {return nullptr;}

protected:
    DataContainer                        M_data;
    SHP(TreeNode)                        M_treeNode;
    SHP(BCManager)                       M_bcManager;
    unsigned int                         M_nComponents;
    EPETRACOMM                           M_comm;
};

}

#include "aAssembler_imp.hpp"

#endif // aASSEMBLER_HPP
