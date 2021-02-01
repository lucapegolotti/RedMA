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
                        shp<TreeNode> treeNode) :
            NavierStokesAssemblerFE(data, treeNode) {
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

        M_TMA_Displacements = TimeMarchingAlgorithmFactory(this->M_data, this->getZeroVector());
        M_TMA_Displacements->setComm(this->M_comm);

        M_boundaryStiffness.reset(new BlockMatrix(this->M_nComponents, this->M_nComponents));
        M_boundaryStiffness->deepCopy(assembleBoundaryStiffness());

        M_boundaryMass.reset(new BlockMatrix(this->M_nComponents, this->M_nComponents));
        M_boundaryMass->deepCopy(assembleBoundaryMass());
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
    assembleBoundaryMass(bool verbose) {
        using namespace LifeV::ExpressionAssembly;

        shp<BlockMatrix> boundaryMass(new BlockMatrix(this->M_nComponents, this->M_nComponents));

        printlog(YELLOW, "Assembling boundary mass matrix ...\n", verbose);

        LifeV::QuadratureBoundary myBDQR(LifeV::buildTetraBDQR(LifeV::quadRuleTria4pt));

        shp<MATRIXEPETRA > BM(new MATRIXEPETRA(M_velocityFESpace->map()));

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

        return boundaryMass;
    }

    shp<aMatrix>
    MembraneAssemblerFE::
    assembleReducedMass(shp<BCManager> bcManager) {
        using namespace LifeV::ExpressionAssembly;

        shp<aMatrix> mass = NavierStokesAssemblerFE::assembleReducedMass(bcManager);

        shp<aMatrix> boundaryMass = this->assembleBoundaryMass();
        boundaryMass->multiplyByScalar(M_membrane_density * M_membrane_thickness +
                                       M_wall_viscoelasticity);

        mass->add(boundaryMass);

        return mass;
    }

    shp<aMatrix>
    MembraneAssemblerFE::
    assembleBoundaryStiffness(bool verbose) {
        using namespace LifeV::ExpressionAssembly;

        shp<BlockMatrix> boundaryStiffness(new BlockMatrix(this->M_nComponents, this->M_nComponents));

        printlog(YELLOW, "Assembling boundary stiffness matrix ...\n", verbose);

        LifeV::QuadratureBoundary myBDQR(LifeV::buildTetraBDQR(LifeV::quadRuleTria4pt));

        shp<MATRIXEPETRA > BS(new MATRIXEPETRA(M_velocityFESpace->map()));

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

        return boundaryStiffness;
    }

    shp<aVector>
    MembraneAssemblerFE::
    getRightHandSide(const double &time,
                     const shp<aVector> &sol_u) {
        shp<aVector> retVec = NavierStokesAssemblerFE::getRightHandSide(time, sol_u);

        shp<aVector> extrapolatedDisplacement = this->M_TMA_Displacements->computeExtrapolatedSolution();

        shp<aVector> membraneContrib = M_boundaryStiffness->multiplyByVector(extrapolatedDisplacement);
        membraneContrib->multiplyByScalar(-M_membrane_thickness);

        shp<aVector> wallContrib = M_boundaryMass->multiplyByVector((extrapolatedDisplacement));
        wallContrib->multiplyByScalar(-M_wall_elasticity);

        retVec->add(membraneContrib);
        retVec->add(wallContrib);

        return retVec;
    }

    void
    MembraneAssemblerFE::
    postProcess(const double &t, const double &dt, const shp<aVector> &sol) {
        NavierStokesAssemblerFE::postProcess(t, dt, sol);

        shp<BlockVector> currDisplacement(new BlockVector(this->M_nComponents));
        currDisplacement->deepCopy(M_TMA_Displacements->simpleAdvance(dt, sol));

        printlog(YELLOW, "Updating displacements field ...\n", this->M_dataContainer.getVerbose());
        M_TMA_Displacements->shiftSolutions(currDisplacement);

        *M_displacementExporter = *static_cast<VECTOREPETRA *>(currDisplacement->block(0)->data().get());
        *M_displacementExporter *= (*M_boundaryIndicator);
    }

    void
    MembraneAssemblerFE::
    setExporter() {
        NavierStokesAssemblerFE::setExporter();

        M_displacementExporter.reset(new VECTOREPETRA(M_velocityFESpace->map(),
                                                      M_exporter->mapType()));

        M_exporter->addVariable(LifeV::ExporterData < MESH > ::VectorField,
                                "displacement", M_velocityFESpace, M_displacementExporter, 0.0);
    }


}  // namespace RedMA
