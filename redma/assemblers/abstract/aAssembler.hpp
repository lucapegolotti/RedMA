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
    aAssembler() {}

    aAssembler(const DataContainer& datafile);

    aAssembler(const DataContainer& datafile, SHP(TreeNode) node);

    inline SHP(BCManager) getBCManager() const {return M_bcManager;}

    virtual void setup() = 0;

    virtual void exportSolution(const double& t,
                                const SHP(aVector)& sol) = 0;

    virtual void postProcess(const double& t,
                             const SHP(aVector)& sol) = 0;

    virtual SHP(aMatrix) getMass(const double& time,
                                const SHP(aVector)& sol) = 0;

    virtual SHP(aMatrix) getMassJacobian(const double& time,
                                         const SHP(aVector)& sol) = 0;

    virtual SHP(aVector) getRightHandSide(const double& time,
                                          const SHP(aVector)& sol) = 0;

    virtual SHP(aMatrix) getJacobianRightHandSide(const double& time,
                                                      const SHP(aVector)& sol) = 0;

    virtual SHP(aVector) getLifting(const double& time) const = 0;

    virtual SHP(aVector) getZeroVector() const = 0;

    virtual void setExporter() = 0;

    virtual void apply0DirichletBCsMatrix(SHP(aMatrix) matrix, double diagCoeff) const = 0;

    virtual void apply0DirichletBCs(SHP(aVector) vector) const = 0;

    virtual void applyDirichletBCs(const double& time, SHP(aVector) vector) const = 0;

    virtual inline SHP(FESPACE) getFEspace(unsigned int index) const {return nullptr;}

    virtual void applyPiola(SHP(aVector) solution, bool inverse) = 0;

    // this must be implemented by the inner assemblers
    virtual inline SHP(FESPACE) getFESpaceBCs() const {return nullptr;}

    virtual inline unsigned int getComponentBCs() const {return 0;}

    virtual inline SHP(TreeNode) getTreeNode() const {return M_treeNode;}

    virtual inline unsigned int getNumComponents() const {return M_nComponents;}

    virtual inline EPETRACOMM getComm() const {return M_comm;}

    virtual inline SHP(ETFESPACE3) getETFESpaceCoupling() const {return nullptr;}

    virtual inline SHP(ETFESPACE1) getETFESpaceSecondary() const {return nullptr;}

    virtual std::vector<SHP(aMatrix)> getMatrices() const {return std::vector<SHP(aMatrix)>();}

    virtual SHP(aMatrix) assembleMatrix(const unsigned int& index,
                                       BlockMDEIMStructure* structure = nullptr) {return SHP(aMatrix)();}

    virtual void setMDEIMs(SHP(MDEIMManager) mdeimManager) {}

    virtual SHP(aVector) getNonLinearTerm() {};

    virtual void initializeFEspaces() {};

    virtual void setDefaultAssemblers(SHP(DefaultAssemblers) defAssemblers)
    {
        M_defaultAssemblers = defAssemblers;
    };

    inline unsigned int ID() {return M_treeNode->M_ID;}

    virtual void RBsetup() {}

    virtual SHP(RBBases) getRBBases() const {}

    virtual void setRBBases(SHP(RBBasesManager) rbManager) {}

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
