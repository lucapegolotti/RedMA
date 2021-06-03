//
// Created by micol on 20/04/2021.
//

#ifndef REDMA_FSIASSEMBLERRB_HPP
#define REDMA_FSIASSEMBLERRB_HPP

#include <redma/assemblers/abstract/aAssemblerFE.hpp>
#include <redma/assemblers/finite_element/FSIAssemblerFE.hpp>
#include <redma/assemblers/reduced_basis/NavierStokesAssemblerRB.hpp>
#include <redma/assemblers/abstract/aAssemblerRB.hpp>
#include <redma/reduced_basis/RBBases.hpp>

namespace RedMA {
    class FSIAssemblerRB : public NavierStokesAssemblerRB {
    public:
        FSIAssemblerRB(const DataContainer &data,
                       shp <TreeNode> treeNode);

        shp <aVector> getRightHandSide(const double &time,
                                       const shp <aVector> &sol) override;

        shp <aMatrix> getJacobianRightHandSide(const double &time,
                                               const shp <aVector> &sol) override;

        void addForcingTermreduced(const shp <aVector> &rhs) const;


        void addFSIMassMatrix(shp<aMatrix> mat);
        virtual void RBsetup() override;

        void postProcess(const double &time, const shp <aVector> &sol) override;
        void exportSolution(const double& t,const shp<aVector>& sol) override;

    protected:
        shp<BlockMatrix>             M_BoundaryStiffnessreduced;


        shp<aTimeMarchingAlgorithm>  M_TMAlgorithm;
        shp<FSIAssemblerFE>          M_FSIAssembler;
    };

}
#endif //REDMA_FSIASSEMBLERRB_HPP
