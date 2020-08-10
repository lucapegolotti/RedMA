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
#include <redma/array/BlockMatrix.hpp>
#include <redma/boundary_conditions/BCManager.hpp>
#include <redma/solver/time_marching_algorithms/aFunctionProvider.hpp>
#include <redma/geometry/TreeStructure.hpp>

#include <redma/reduced_basis/MDEIMStructure.hpp>
#include <redma/reduced_basis/MDEIMManager.hpp>
#include <redma/reduced_basis/RBBasesManager.hpp>

#include <redma/problem/DataContainer.hpp>
#include <redma/assemblers/DefaultAssemblersLibrary.hpp>

namespace RedMA
{

class aAssembler : public aFunctionProvider
{
    typedef DefaultAssemblersLibrary DefaultAssemblers;
public:
    aAssembler(const DataContainer& datafile);

    aAssembler(const DataContainer& datafile, SHP(TreeNode) node);

    inline SHP(BCManager) getBCManager() const {return M_bcManager;}

    virtual void setup() = 0;

    virtual void exportSolution(const double& t,
                                const BlockVector& sol) = 0;

    virtual void postProcess(const double& t,
                             const BlockVector& sol) = 0;

    virtual BlockMatrix getMass(const double& time,
                                const BlockVector& sol) = 0;

    virtual BlockMatrix getMassJacobian(const double& time,
                                        const BlockVector& sol) = 0;

    virtual BlockVector getRightHandSide(const double& time,
                                         const BlockVector& sol) = 0;

    virtual BlockMatrix getJacobianRightHandSide(const double& time,
                                                 const BlockVector& sol) = 0;

    virtual BlockVector getLifting(const double& time) const = 0;

    virtual BlockVector getFELifting(const double& time) const {}

    virtual BlockVector getZeroVector() const = 0;

    virtual void setExporter() = 0;

    virtual void apply0DirichletBCsMatrix(BlockMatrix& matrix, double diagCoeff) const = 0;

    virtual void apply0DirichletBCs(BlockVector& vector) const = 0;

    virtual void applyDirichletBCs(const double& time, BlockVector& vector) const = 0;

    virtual inline SHP(FESPACE) getFEspace(unsigned int index) const {return nullptr;}

    virtual void applyPiola(BlockVector solution, bool inverse) = 0;

    // this must be implemented by the inner assemblers
    virtual inline SHP(FESPACE) getFESpaceBCs() const {return nullptr;}

    virtual inline unsigned int getComponentBCs() const {return 0;}

    virtual inline SHP(TreeNode) getTreeNode() const {return M_treeNode;}

    virtual inline unsigned int getNumComponents() const {return M_nComponents;}

    virtual inline EPETRACOMM getComm() const {return M_comm;}

    virtual inline SHP(ETFESPACE3) getETFESpaceCoupling() const {return nullptr;}

    virtual inline SHP(ETFESPACE1) getETFESpaceSecondary() const {return nullptr;}

    virtual std::vector<BlockMatrix> getMatrices() const {return std::vector<BlockMatrix>();}

    virtual BlockMatrix assembleMatrix(const unsigned int& index,
                                       BlockMDEIMStructure* structure = nullptr) {return BlockMatrix();}

    virtual void setMDEIMs(SHP(MDEIMManager) mdeimManager) {}

    virtual BlockVector getNonLinearTerm() {};

    virtual void initializeFEspaces() {};

    virtual void setDefaultAssemblers(SHP(DefaultAssemblers) defAssemblers)
    {
        M_defaultAssemblers = defAssemblers;
    };

    inline unsigned int ID() {return M_treeNode->M_ID;}

protected:
    DataContainer                           M_data;
    SHP(TreeNode)                           M_treeNode;
    SHP(BCManager)                          M_bcManager;
    unsigned int                            M_nComponents;
    EPETRACOMM                              M_comm;
    std::string                             M_name;
    SHP(DefaultAssemblers)                  M_defaultAssemblers;
};

}

#endif // aASSEMBLER_HPP
