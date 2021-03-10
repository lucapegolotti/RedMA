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
        M_transverse_shear_coeff = data("structure/transverse_shear_coefficient", 1.0);
        M_wall_elasticity = data("structure/external_wall/elastic", 1e4);
        M_wall_viscoelasticity = data("structure/external_wall/plastic", 1e3);
        M_wallFlag = data("structure/flag", 10);

        this->computeLameConstants();

        M_addNoSlipBC = false;
    }

    void
    MembraneAssemblerFE::
    setup() {

        // computing thickness in advance, as it is needed to assemble boundary matrices
        this->computeThickness();

        NavierStokesAssemblerFE::setup();

        std::vector<unsigned int> flags = M_bcManager->getWallFlags(true);
        M_boundaryIndicator = M_bcManager->computeBoundaryIndicator(M_velocityFESpace, flags);

        if (!(M_TMA_Displacements)){
            M_TMA_Displacements = TimeMarchingAlgorithmFactory(this->M_data, this->getZeroVector());
            M_TMA_Displacements->setComm(this->M_comm);
        }

        // exporting thickness as scalar field at time t0
        this->exportThickness();
    }

    void
    MembraneAssemblerFE::
    computeLameConstants() {
        double poisson = this->M_data("structure/poisson", 0.30);
        double young = this->M_data("structure/young", 4e6);

        // this is not exactly Lame1 (i.e. lambda), but it is the coefficient appearing in the membrane model
        M_lameI = (young * poisson) / ((1. - poisson) * (1. + poisson));
        M_lameII = young / (2. * (1. + poisson));

    }

    void
    MembraneAssemblerFE::
    computeThickness()
    {
        this->M_treeNode->M_block->computeMembraneThickness();
    }

    void
    MembraneAssemblerFE::
    exportThickness()
    {
        if (!M_exporter)
            StokesAssemblerFE::setExporter();

        double t0 = M_data("time_discretization/t0", 0.0);

        CoutRedirecter ct1;
        ct1.redirect();

        // setting thickness to 0 in the interior of the domain
        std::vector<unsigned int> flags = M_bcManager->getWallFlags(true);
        shp<VECTOREPETRA> boundaryIndicator = M_bcManager->computeBoundaryIndicator(M_pressureFESpace,
                                                                                  flags);
        shp<VECTOREPETRA> thickness = this->M_treeNode->M_block->getMembraneThickness();
        *thickness *= (*boundaryIndicator);

        M_exporter->addVariable(LifeV::ExporterData<MESH>::ScalarField,
                                "thickness", M_pressureFESpace,
                                thickness, 0.0);

        printlog(CYAN, ct1.restore(), false);

        CoutRedirecter ct2;
        ct2.redirect();

        M_exporter->postProcess(t0);

        printlog(CYAN, ct2.restore());
    }

    shp<aMatrix>
    MembraneAssemblerFE::
    assembleBoundaryMass(shp<BCManager> bcManager, bool verbose) {
        using namespace LifeV::ExpressionAssembly;

        shp<BlockMatrix> boundaryMass(new BlockMatrix(this->M_nComponents, this->M_nComponents));

        printlog(YELLOW, "Assembling boundary mass matrix ...\n", verbose);

        LifeV::QuadratureBoundary myBDQR(LifeV::buildTetraBDQR(LifeV::quadRuleTria4pt));

        shp<VECTOREPETRA> thickness = this->M_treeNode->M_block->getMembraneThickness();
        shp<VECTOREPETRA>  thicknessRepeated(new VECTOREPETRA(*thickness, LifeV::Repeated));

        shp<MATRIXEPETRA> BM(new MATRIXEPETRA(M_velocityFESpace->map()));
        integrate(boundary(M_velocityFESpaceETA->mesh(), M_wallFlag),
                  myBDQR,
                  M_velocityFESpaceETA,
                  M_velocityFESpaceETA,
                  value(M_membrane_density) *
                  value(M_pressureFESpaceETA, *thicknessRepeated) *
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
    assembleWallBoundaryMass(shp<BCManager> bcManager, bool verbose) {
        using namespace LifeV::ExpressionAssembly;

        shp<BlockMatrix> wallBoundaryMass(new BlockMatrix(this->M_nComponents, this->M_nComponents));

        printlog(YELLOW, "Assembling boundary mass matrix ...\n", verbose);

        LifeV::QuadratureBoundary myBDQR(LifeV::buildTetraBDQR(LifeV::quadRuleTria4pt));

        shp<MATRIXEPETRA> WBM(new MATRIXEPETRA(M_velocityFESpace->map()));

        integrate(boundary(M_velocityFESpaceETA->mesh(), M_wallFlag),
                  myBDQR,
                  M_velocityFESpaceETA,
                  M_velocityFESpaceETA,
                  dot(phi_i, phi_j)
        ) >> WBM;
        WBM->globalAssemble();

        shp<SparseMatrix> Mwrapper(new SparseMatrix);
        Mwrapper->setData(WBM);
        wallBoundaryMass->setBlock(0, 0, Mwrapper);

        bcManager->apply0DirichletMatrix(*wallBoundaryMass, M_velocityFESpace,
                                         0, 0.0, !(M_addNoSlipBC));

        return wallBoundaryMass;
    }

    shp<aMatrix>
    MembraneAssemblerFE::
    assembleMass(shp<BCManager> bcManager) {
        using namespace LifeV::ExpressionAssembly;

        shp<aMatrix> mass = NavierStokesAssemblerFE::assembleMass(bcManager);

        M_boundaryMass.reset(new BlockMatrix(this->M_nComponents, this->M_nComponents));
        M_boundaryMass->deepCopy(this->assembleBoundaryMass(M_bcManager));

        M_wallBoundaryMass.reset(new BlockMatrix(this->M_nComponents, this->M_nComponents));
        M_wallBoundaryMass->deepCopy(assembleWallBoundaryMass(M_bcManager));

        mass->add(M_boundaryMass);

        return mass;
    }

    shp<aMatrix>
    MembraneAssemblerFE::
    assembleBoundaryStiffness(shp<BCManager> bcManager, bool verbose) {
        using namespace LifeV::ExpressionAssembly;

        shp<BlockMatrix> boundaryStiffness(new BlockMatrix(this->M_nComponents, this->M_nComponents));

        printlog(YELLOW, "Assembling boundary stiffness matrix ...\n", verbose);

        LifeV::QuadratureBoundary myBDQR(LifeV::buildTetraBDQR(LifeV::quadRuleTria4pt));

        shp<VECTOREPETRA> thickness = this->M_treeNode->M_block->getMembraneThickness();
        shp<VECTOREPETRA>  thicknessRepeated(new VECTOREPETRA(*thickness, LifeV::Repeated));

        LifeV::MatrixSmall<3, 3> Eye;
        Eye *= 0.0;
        Eye[0][0] = 1;
        Eye[1][1] = 1;
        Eye[2][2] = 1;

        shp<MATRIXEPETRA> BS(new MATRIXEPETRA(M_velocityFESpace->map()));
        integrate(boundary(M_velocityFESpaceETA->mesh(), M_wallFlag),
                  myBDQR,
                  M_velocityFESpaceETA,
                  M_velocityFESpaceETA,
                  value(M_pressureFESpaceETA, *thicknessRepeated) * (
                  2.0 * value(this->M_lameII) *
                  0.5 * dot(
                          (grad(phi_j) - grad(phi_j) * outerProduct(Nface, Nface))
                          + transpose(grad(phi_j) - grad(phi_j) * outerProduct(Nface, Nface)),
                          (grad(phi_i) - grad(phi_i) * outerProduct(Nface, Nface))) +
                  value(this->M_lameI) *
                          dot(value(Eye), (grad(phi_j) - grad(phi_j) * outerProduct(Nface, Nface))) *
                          dot(value(Eye), (grad(phi_i) - grad(phi_i) * outerProduct(Nface, Nface))) +
                  2.0 * value((this->M_transverse_shear_coeff - 1.0) * this->M_lameII) *
                  0.5 * dot(
                          ((grad(phi_j) - grad(phi_j) * outerProduct(Nface, Nface) +
                          transpose(grad(phi_j) - grad(phi_j) * outerProduct(Nface, Nface))) *
                          outerProduct(Nface, Nface)),
                          (grad(phi_i) - grad(phi_i) * outerProduct(Nface, Nface)))
                                                                    )
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
    assembleStiffness(shp<BCManager> bcManager) {
        using namespace LifeV::ExpressionAssembly;

        shp<aMatrix> stiffness = NavierStokesAssemblerFE::assembleStiffness(bcManager);

        if (!(M_TMA_Displacements)){
            M_TMA_Displacements = TimeMarchingAlgorithmFactory(this->M_data, this->getZeroVector());
            M_TMA_Displacements->setComm(this->M_comm);
        }

        double  dt = this->M_data("time_discretization/dt", 0.01);
        double rhs_coeff = this->M_TMA_Displacements->getCoefficients().back();

        M_boundaryStiffness.reset(new BlockMatrix(this->M_nComponents, this->M_nComponents));
        M_boundaryStiffness->deepCopy(assembleBoundaryStiffness(M_bcManager));
        M_boundaryStiffness->multiplyByScalar(dt * rhs_coeff);

        stiffness->add(M_boundaryStiffness);

        M_boundaryStiffness->multiplyByScalar(1.0 / (dt * rhs_coeff));

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

        // computing external wall contribution involved as velocity mass matrix in system matrix
        double  dt = this->M_data("time_discretization/dt", 0.01);
        double rhs_coeff = this->M_TMA_Displacements->getCoefficients().back();

        shp<BlockMatrix> wallMatrix(new BlockMatrix(this->M_nComponents,
                                                      this->M_nComponents));

        wallMatrix->add(this->M_wallBoundaryMass);
        wallMatrix->multiplyByScalar(-1.0 * M_wall_viscoelasticity
                                     - 1.0 * M_wall_elasticity * dt * rhs_coeff);

        // computing current displacement part due to previous displacements
        shp<aVector> rhsDisplacement = this->M_TMA_Displacements->combineOldSolutions();

        // computing membrane stress contribution on previous displacements
        shp<aVector> membraneContrib = M_boundaryStiffness->multiplyByVector(rhsDisplacement);

        // computing external wall contribution on previous displacements
        shp<aVector> wallContrib = M_wallBoundaryMass->multiplyByVector((rhsDisplacement));
        wallContrib->multiplyByScalar(M_wall_elasticity);

        // adding the three new contributions
        retVec->add(wallMatrix->multiplyByVector(sol));
        retVec->add(membraneContrib);
        retVec->add(wallContrib);

        this->M_bcManager->apply0DirichletBCs(*spcast<BlockVector>(retVec), this->getFESpaceBCs(),
                                              this->getComponentBCs(), !(this->M_addNoSlipBC));

        return retVec;
    }

    shp<aMatrix>
    MembraneAssemblerFE::
    getJacobianRightHandSide(const double &time, const shp<aVector> &sol)
    {
        shp<aMatrix> retMat = NavierStokesAssemblerFE::getJacobianRightHandSide(time, sol);

        // adding external wall contribution involved as velocity mass matrix in system matrix
        double  dt = this->M_data("time_discretization/dt", 0.01);
        double rhs_coeff = this->M_TMA_Displacements->getCoefficients().back();

        shp<BlockMatrix> wallMatrix(new BlockMatrix(this->M_nComponents,
                                                      this->M_nComponents));

        wallMatrix->add(this->M_wallBoundaryMass);
        wallMatrix->multiplyByScalar(-1.0 * M_wall_viscoelasticity
                                     - 1.0 * M_wall_elasticity * dt * rhs_coeff);

        retMat->add(wallMatrix);

        this->M_bcManager->apply0DirichletMatrix(*spcast<BlockMatrix>(retMat),
                                                 this->getFESpaceBCs(),
                                                 this->getComponentBCs(), 0.0,
                                                 !(this->M_addNoSlipBC));

        return retMat;
    }

    void
    MembraneAssemblerFE::
    postProcess(const double &t, const double &dt, const shp<aVector> &sol) {
        NavierStokesAssemblerFE::postProcess(t, dt, sol);

        printlog(YELLOW, "[MembraneAssemblerFE] Updating displacements field ...\n",
                 this->M_data.getVerbose());

        shp<BlockVector> currDisplacement(new BlockVector(this->M_nComponents));
        currDisplacement->deepCopy(M_TMA_Displacements->simpleAdvance(dt, convert<BlockVector>(sol)));

        // TODO: it may be a good idea to export in exportSolution, but I should well figure out how
        // TODO: to do it in a smart way!
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

        M_exporter->addVariable(LifeV::ExporterData<MESH>::VectorField,
                                "displacement", M_velocityFESpace, M_displacementExporter, 0.0);
    }

}  // namespace RedMA
