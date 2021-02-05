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

#include "MembraneAssemblerFE.hpp"

namespace RedMA {

    MembraneAssemblerFE::
    MembraneAssemblerFE(const DataContainer &data,
                        shp<TreeNode> treeNode,
                        std::string stabilizationName) :
            NavierStokesAssemblerFE(data, treeNode, stabilizationName) {
        M_membrane_density = data("structure/density", 1.2);
        M_membrane_thickness = data("structure/thickness", 0.1);
        M_wall_elasticity = data("structure/external_wall/elastic", 1e4);
        M_wall_viscoelasticity = data("structure/external_wall/plastic", 1e3);
        M_wallFlag = data("structure/flag", 10);

        computeLameConstants();

        M_addNoSlipBC = false;
    }

    void
    MembraneAssemblerFE::
    setup() {

        NavierStokesAssemblerFE::setup();

        M_boundaryIndicator = M_bcManager->computeBoundaryIndicator(M_velocityFESpace);

        if (!(M_TMA_Displacements)){
            M_TMA_Displacements = TimeMarchingAlgorithmFactory(this->M_data, this->getZeroVector());
            M_TMA_Displacements->setComm(this->M_comm);
        }

        M_boundaryStiffness.reset(new BlockMatrix(this->M_nComponents, this->M_nComponents));
        M_boundaryStiffness->deepCopy(assembleBoundaryStiffness(M_bcManager));

        M_boundaryMass.reset(new BlockMatrix(this->M_nComponents, this->M_nComponents));
        M_boundaryMass->deepCopy(assembleBoundaryMass(M_bcManager));
    }

    void
    MembraneAssemblerFE::
    computeLameConstants() {
        double poisson = this->M_dataContainer("structure/poisson", 0.45);
        double young = this->M_dataContainer("structure/young", 4e6);

        M_lameI = (young * poisson) / ((1. - 2 * poisson) * (1. + poisson));
        M_lameII = young / (2. * (1. + poisson));

    }

    shp<aMatrix>
    MembraneAssemblerFE::
    assembleBoundaryMass(shp<BCManager> bcManager, bool verbose) {
        using namespace LifeV::ExpressionAssembly;

        shp<BlockMatrix> boundaryMass(new BlockMatrix(this->M_nComponents, this->M_nComponents));

        printlog(YELLOW, "Assembling boundary mass matrix ...\n", verbose);

        LifeV::QuadratureBoundary myBDQR(LifeV::buildTetraBDQR(LifeV::quadRuleTria4pt));

        shp<MATRIXEPETRA> BM(new MATRIXEPETRA(M_velocityFESpace->map()));

        integrate(boundary(M_velocityFESpaceETA->mesh(), M_wallFlag),
                  myBDQR,
                  M_velocityFESpaceETA,
                  M_velocityFESpaceETA,
                  dot(phi_i, phi_j)
        ) >> BM;
        BM->globalAssemble();

        shp<SparseMatrix> Mwrapper(new SparseMatrix);
        Mwrapper->setData(BM);
        boundaryMass->setBlock(0, 0, Mwrapper);

        bcManager->apply0DirichletMatrix(*boundaryMass, M_velocityFESpace, 0, 0.0, !(M_addNoSlipBC));

        return boundaryMass;
    }

    shp<aMatrix>
    MembraneAssemblerFE::
    assembleReducedMass(shp<BCManager> bcManager) {
        using namespace LifeV::ExpressionAssembly;

        shp<aMatrix> mass = NavierStokesAssemblerFE::assembleReducedMass(bcManager);

        shp<aMatrix> boundaryMass = this->assembleBoundaryMass(bcManager);
        boundaryMass->multiplyByScalar(M_membrane_density * M_membrane_thickness);

        mass->add(boundaryMass);

        return mass;
    }

    shp<aMatrix>
    MembraneAssemblerFE::
    assembleBoundaryStiffness(shp<BCManager> bcManager, bool verbose) {
        using namespace LifeV::ExpressionAssembly;

        shp<BlockMatrix> boundaryStiffness(new BlockMatrix(this->M_nComponents, this->M_nComponents));

        printlog(YELLOW, "Assembling boundary stiffness matrix ...\n", verbose);

        LifeV::QuadratureBoundary myBDQR(LifeV::buildTetraBDQR(LifeV::quadRuleTria4pt));

        shp<MATRIXEPETRA> BS(new MATRIXEPETRA(M_velocityFESpace->map()));

        LifeV::MatrixSmall<3, 3> Eye;
        Eye *= 0.0;
        Eye[0][0] = 1;
        Eye[1][1] = 1;
        Eye[2][2] = 1;

        integrate(boundary(M_velocityFESpaceETA->mesh(), M_wallFlag),
                  myBDQR,
                  M_velocityFESpaceETA,
                  M_velocityFESpaceETA,
                  2 * value(this->M_lameII) *
                  0.5 * dot(
                          (grad(phi_j) - grad(phi_j) * outerProduct(Nface, Nface))
                          + transpose(grad(phi_j) - grad(phi_j) * outerProduct(Nface, Nface)),
                          (grad(phi_i) - grad(phi_i) * outerProduct(Nface, Nface))) +
                  value(this->M_lameI) *
                  dot(value(Eye), (grad(phi_j) - grad(phi_j) * outerProduct(Nface, Nface))) *
                  dot(value(Eye), (grad(phi_i) - grad(phi_i) * outerProduct(Nface, Nface)))
        ) >> BS;
        BS->globalAssemble();

        shp<SparseMatrix> Awrapper(new SparseMatrix);
        Awrapper->setData(BS);
        boundaryStiffness->setBlock(0, 0, Awrapper);

        bcManager->apply0DirichletMatrix(*boundaryStiffness, M_velocityFESpace, 0, 0.0, !(M_addNoSlipBC));

        return boundaryStiffness;
    }

    shp<aMatrix>
    MembraneAssemblerFE::
    assembleReducedStiffness(shp<BCManager> bcManager) {
        using namespace LifeV::ExpressionAssembly;

        shp<aMatrix> stiffness = NavierStokesAssemblerFE::assembleReducedStiffness(bcManager);

        if (!(M_TMA_Displacements)){
            M_TMA_Displacements = TimeMarchingAlgorithmFactory(this->M_data, this->getZeroVector());
            M_TMA_Displacements->setComm(this->M_comm);
        }

        double  dt = this->M_dataContainer("time_discretization/dt", 0.01);
        double rhs_coeff = this->M_TMA_Displacements->getCoefficients().back();

        shp<aMatrix> boundaryStiffness = this->assembleBoundaryStiffness(bcManager);
        boundaryStiffness->multiplyByScalar(M_membrane_thickness * dt * rhs_coeff);

        stiffness->add(boundaryStiffness);

        return stiffness;
    }

    shp<aVector>
    MembraneAssemblerFE::
    getRightHandSide(const double &time,
                     const shp<aVector> &sol) {
        shp<aVector> retVec = NavierStokesAssemblerFE::getRightHandSide(time, sol);

        if (!(M_TMA_Displacements)){
            M_TMA_Displacements = TimeMarchingAlgorithmFactory(this->M_data, this->getZeroVector());
            M_TMA_Displacements->setComm(this->M_comm);
        }

        // adding external wall contribution involved as velocity mass matrix in system matrix
        double  dt = this->M_dataContainer("time_discretization/dt", 0.01);
        double rhs_coeff = this->M_TMA_Displacements->getCoefficients().back();

        shp<BlockMatrix> wallMatrix(new BlockMatrix(this->M_nComponents,
                                                      this->M_nComponents));

        wallMatrix->add(M_boundaryMass);
        wallMatrix->multiplyByScalar(M_wall_viscoelasticity + M_wall_elasticity * dt * rhs_coeff);

        retVec->add(wallMatrix->multiplyByVector(sol));

        // computing current displacement part due to previous displacements
        shp<aVector> rhsDisplacement = this->M_TMA_Displacements->combineOldSolutions();

        // computing membrane stress contribution on previous displacements
        shp<aVector> membraneContrib = M_boundaryStiffness->multiplyByVector(rhsDisplacement);
        membraneContrib->multiplyByScalar(M_membrane_thickness);

        // computing external wall contribution on previous displacements
        shp<aVector> wallContrib = M_boundaryMass->multiplyByVector((rhsDisplacement));
        wallContrib->multiplyByScalar(M_wall_elasticity);

        retVec->add(membraneContrib);
        retVec->add(wallContrib);

        this->M_bcManager->apply0DirichletBCs(*spcast<BlockVector>(retVec), this->getFESpaceBCs(),
                                              this->getComponentBCs(), !(this->M_addNoSlipBC));

        return retVec;
    }

    void
    MembraneAssemblerFE::
    postProcess(const double &t, const double &dt, const shp<aVector> &sol) {
        NavierStokesAssemblerFE::postProcess(t, dt, sol);

        printlog(YELLOW, "[MembraneAssemblerFE] Updating displacements field ...\n",
                 this->M_dataContainer.getVerbose());

        shp<BlockVector> currDisplacement(new BlockVector(this->M_nComponents));
        currDisplacement->deepCopy(M_TMA_Displacements->simpleAdvance(dt, convert<BlockVector>(sol)));

        *M_displacementExporter = *static_cast<VECTOREPETRA *>(currDisplacement->block(0)->data().get());
        *M_displacementExporter *= (*M_boundaryIndicator);

        M_TMA_Displacements->shiftSolutions(currDisplacement);
    }

    void
    MembraneAssemblerFE::
    setExporter() {
        NavierStokesAssemblerFE::setExporter();

        M_displacementExporter.reset(new VECTOREPETRA(M_velocityFESpace->map(),
                                                      M_exporter->mapType()));

        M_exporter->addVariable(LifeV::ExporterData <MESH> ::VectorField,
                                "displacement", M_velocityFESpace, M_displacementExporter, 0.0);
    }


}  // namespace RedMA
