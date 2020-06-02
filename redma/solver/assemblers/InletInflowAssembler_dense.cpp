#include "InletInflowAssembler.hpp"

namespace RedMA
{

template <>
void
InletInflowAssembler<DenseVector, DenseMatrix>::
addContributionRhs(const double& time,
                   BlockVector<BlockVector<DenseVector>>& rhs,
                   const BlockVector<BlockVector<DenseVector>>& sol,
                   const unsigned int& nPrimalBlocks)
{
    unsigned int childID = this->M_interface.M_indexChild;
    unsigned int interfaceID = this->M_interface.M_ID;
    SHP(aAssembler<DenseVector COMMA DenseMatrix>) assemblerChild;
    assemblerChild = this->M_interface.M_assemblerChild;

    rhs.block(childID) -= this->M_childBT * sol.block(nPrimalBlocks + interfaceID);

    rhs.block(nPrimalBlocks + interfaceID) -= this->M_childB * sol.block(childID);
    BlockVector<VectorEp> lifting = assemblerChild->getFELifting(time);
    auto fecontribution = this->M_childBEp * lifting;
    auto projection = assemblerChild->getRBBases()->projectOnLagrangeSpace(fecontribution);
    rhs.block(nPrimalBlocks + interfaceID) += projection;
    // BlockVector<VectorEp> lifting = assemblerChild->getFELifting(time);
    // BlockVector<DenseVector> liftingProjected = assemblerChild->getRBBases()->leftProject(lifting);
    // rhs.block(nPrimalBlocks + interfaceID) += this->M_childB * liftingProjected;
}
//
// template <class InVectorType, class InMatrixType>
// void
// InletInflowAssembler<InVectorType, InMatrixType>::
// addContributionJacobianRhs(const double& time,
//                            BlockMatrix<BlockMatrix<InMatrixType>>& jac,
//                            const BlockVector<BlockVector<InVectorType>>& sol,
//                            const unsigned int& nPrimalBlocks)
// {
//     unsigned int fatherID = this->M_interface.M_indexFather;
//     unsigned int childID = this->M_interface.M_indexChild;
//     unsigned int interfaceID = this->M_interface.M_ID;
//
//     // hard copy, otherwise we flip the sign of the matrices every time this
//     // function is called
//     jac.block(childID,  nPrimalBlocks + interfaceID).hardCopy(this->M_childBT);
//     jac.block(nPrimalBlocks + interfaceID,  childID).hardCopy(this->M_childB);
//
//     jac.block(childID,  nPrimalBlocks + interfaceID) *= (-1);
//     jac.block(nPrimalBlocks + interfaceID,  childID) *= (-1);
//
//     if (this->M_stabilizationCoupling > THRESHOLDSTAB)
//     {
//         jac.block(nPrimalBlocks + interfaceID,  childID) += (this->M_stabChild * (-1.0 * this->M_stabilizationCoupling));
//
//         jac.block(nPrimalBlocks + interfaceID, nPrimalBlocks + interfaceID).hardCopy(this->M_identity * (-1.0 * this->M_stabilizationCoupling));
//     }
// }

}
