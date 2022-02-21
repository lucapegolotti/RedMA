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

#include "MembraneAssemblerRB.hpp"

namespace RedMA {

MembraneAssemblerRB::
MembraneAssemblerRB(const DataContainer& data, shp<TreeNode> treeNode):
    NavierStokesAssemblerRB(data, treeNode)
{
    replaceFEAssembler(std::make_shared<MembraneAssemblerFE>(data, treeNode));
    M_name = "MembraneAssemblerRB";
}

void
MembraneAssemblerRB::
RBsetup()
{
    NavierStokesAssemblerRB::RBsetup();

    unsigned int id = M_treeNode->M_ID;

    printlog(YELLOW, "[MembraneAssemblerRB] assembling and projecting boundary matrices",
             this->M_data.getVerbose());
    Chrono chrono;
    chrono.start();

    M_reducedBoundaryMass.reset(new BlockMatrix(2,2));
    M_reducedBoundaryStiffness.reset(new BlockMatrix(2,2));
    M_reducedWallBoundaryMass.reset(new BlockMatrix(2,2));

    auto matrices = this->getFEAssembler()->getBoundaryMatrices();
    M_reducedBoundaryMass->setBlock(0,0,M_bases->matrixProject(matrices[0]->block(0,0),
                                                                   0, 0, id));
    M_reducedBoundaryStiffness->setBlock(0,0,M_bases->matrixProject(matrices[1]->block(0,0),
                                                                   0, 0, id));
    M_reducedWallBoundaryMass->setBlock(0,0,M_bases->matrixProject(matrices[2]->block(0,0),
                                                                   0, 0, id));

    std::string msg = "done, in ";
    msg += std::to_string(chrono.diff());
    msg += " seconds\n";
    printlog(YELLOW, msg, this->M_data.getVerbose());

    if (!(M_TMA_Displacements)){
        M_TMA_Displacements = TimeMarchingAlgorithmFactory(this->M_data, this->getZeroVector());
        M_TMA_Displacements->setComm(this->M_comm);
    }

}

shp<aVector>
MembraneAssemblerRB::
getRightHandSide(const double& time,
                 const shp<aVector>& sol)
{
    Chrono chrono;
    chrono.start();

    shp<BlockVector> solBlck = convert<BlockVector>(sol);

    std::string msg = "[MembraneAssemblerRB] computing right-hand side term ...";
    printlog(YELLOW, msg, this->M_data.getVerbose());

    shp<BlockVector> retVec = convert<BlockVector>(NavierStokesAssemblerRB::getRightHandSide(time,sol));

    if (!(M_TMA_Displacements)){
        M_TMA_Displacements = TimeMarchingAlgorithmFactory(this->M_data, this->getZeroVector());
        M_TMA_Displacements->setComm(this->M_comm);
    }

    // computing external wall contribution involved as velocity mass matrix in system matrix
    double  dt = this->M_data("time_discretization/dt", 0.01);
    // TODO: here we can find a way of not making any interaction with the time marching scheme!
    double rhs_coeff = this->M_TMA_Displacements->getCoefficients().back();

    shp<BlockMatrix> wallMatrix(new BlockMatrix(this->M_nComponents,
                                                  this->M_nComponents));

    wallMatrix->add(this->M_reducedWallBoundaryMass);
    wallMatrix->multiplyByScalar(-1.0 * this->getFEAssembler()->getWallViscoelasticity()
                                 - 1.0 * this->getFEAssembler()->getWallElasticity() * dt * rhs_coeff);

    // computing current displacement part due to previous displacements
    shp<aVector> rhsDisplacement = this->M_TMA_Displacements->combineOldSolutions();

    // computing membrane stress contribution on previous displacements
    shp<aVector> membraneContrib = M_reducedBoundaryStiffness->multiplyByVector(rhsDisplacement);

    // computing external wall contribution on previous displacements
    shp<aVector> wallContrib = M_reducedWallBoundaryMass->multiplyByVector((rhsDisplacement));
    wallContrib->multiplyByScalar(this->getFEAssembler()->getWallElasticity());

    // adding the three new contributions
    retVec->add(wallMatrix->multiplyByVector(sol));
    retVec->add(membraneContrib);
    retVec->add(wallContrib);

    msg = "done, in ";
    msg += std::to_string(chrono.diff());
    msg += " seconds\n";
    printlog(YELLOW, msg, this->M_data.getVerbose());

    return retVec;
}

shp<aMatrix>
MembraneAssemblerRB::
getJacobianRightHandSide(const double &time,
                         const shp<aVector> &sol,
                         const double& diagCoeff)
{
    Chrono chrono;
    chrono.start();

    std::string msg = "[MembraneAssemblerRB] computing right-hand side jacobian ...";
    printlog(YELLOW, msg, this->M_data.getVerbose());

    shp<aMatrix> retMat = NavierStokesAssemblerRB::getJacobianRightHandSide(time, sol);

    // adding external wall contribution involved as velocity mass matrix in system matrix
    double  dt = this->M_data("time_discretization/dt", 0.01);
    // TODO: here we can find a way of not making any interaction with the time marching scheme!
    double rhs_coeff = this->M_TMA_Displacements->getCoefficients().back();

    shp<BlockMatrix> wallMatrix(new BlockMatrix(this->M_nComponents,
                                                  this->M_nComponents));

    wallMatrix->add(this->M_reducedWallBoundaryMass);
    wallMatrix->multiplyByScalar(-1.0 * this->getFEAssembler()->getWallViscoelasticity()
                                 - 1.0 * this->getFEAssembler()->getWallElasticity() * dt * rhs_coeff);

    retMat->add(wallMatrix);

    msg = "done, in ";
    msg += std::to_string(chrono.diff());
    msg += " seconds\n";
    printlog(YELLOW, msg, this->M_data.getVerbose());

    return retMat;
}

void
MembraneAssemblerRB::
postProcess(const double &t, const shp<aVector> &sol)
{
    StokesAssemblerRB::postProcess(t, sol);

    printlog(YELLOW, "[MembraneAssemblerRB] Updating displacements field ...\n",
             this->M_data.getVerbose());

    unsigned int id = M_treeNode->M_ID;
    double  dt = this->M_data("time_discretization/dt", 0.01);

    shp<BlockVector> currDisplacement(new BlockVector(this->M_nComponents));
    currDisplacement->deepCopy(M_TMA_Displacements->simpleAdvance(dt, convert<BlockVector>(sol)));

    // TODO: it may be a good idea to export in exportSolution, but I should well figure out how
    // TODO: to do it in a smart way!
    shp<VECTOREPETRA> exportedDisplacement = M_bases->reconstructFEFunction(currDisplacement->block(0), 0, id);
    *exportedDisplacement *= *(this->getFEAssembler()->getBoundaryIndicator());
    this->getFEAssembler()->setDisplacementExporter(exportedDisplacement);

    M_TMA_Displacements->shiftSolutions(currDisplacement);
}

}  // namespace RedMA
