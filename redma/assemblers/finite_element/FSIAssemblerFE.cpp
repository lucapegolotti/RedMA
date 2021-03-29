//
// Created by micol on 13/03/2021.
//
#include "FSIAssemblerFE.hpp"

namespace RedMA
{
  FSIAssemblerFE::
  FSIAssemblerFE(const DataContainer& data, shp<TreeNode> treeNode,std::string stabilizationName):
                          NavierStokesAssemblerFE(data, treeNode, stabilizationName),M_BoundaryStiffness(new BlockMatrix(this->M_nComponents,this->M_nComponents))
                          ,M_TMAlgorithm(TimeMarchingAlgorithmFactory(data)){


  }
  void
  FSIAssemblerFE::computeLameConstants(){
    double poisson = M_data("structure/poisson", 0.45);
    double young = M_data("structure/young", 4e6);
    double thickness = M_data("structure/thickness", 0.1);
    M_lameI = (thickness * young * poisson)/((1. - poisson)*(1. + poisson));
    M_lameII = thickness * young/(2. * (1. + poisson));
  }

  void
  FSIAssemblerFE::setup(){
      NavierStokesAssemblerFE::setup();
      //initialize BDF with the initial displacement??
      M_TMAlgorithm->setComm(M_comm);
      M_TMAlgorithm->setup(this->getZeroVector());
      this->computeLameConstants();
      this->getBoundaryStiffnessMatrix();

      this->addFSIMassMatrix(M_mass);
      this->M_bcManager->apply0DirichletMatrix(*M_mass,this->getFESpaceBCs(),this->getComponentBCs(), 1.0);

      printlog(YELLOW, "[FSI]Mass matrix assembled ...\n", this->M_data.getVerbose());
      this->setExporterDisplacement();
  }
  void
  FSIAssemblerFE::
  addFSIMassMatrix(shp<aMatrix> mat)
    {
        using namespace LifeV;
        using namespace ExpressionAssembly;

        printlog(YELLOW, "[NS]Assembling boundary mass matrix ...\n", this->M_data.getVerbose());

        LifeV::QuadratureBoundary myBDQR(LifeV::buildTetraBDQR(LifeV::quadRuleTria4pt)); //not sure if updated, there is another way
        shp<MATRIXEPETRA> M_boundaryMass(new MATRIXEPETRA(M_velocityFESpace->map()));

        double density = M_data("structure/density", 1.2);
        double thickness = M_data("structure/thickness", 0.1);
        unsigned int wallFlag = M_data("structure/flag", 10);
        integrate(boundary(M_velocityFESpaceETA->mesh(), wallFlag),
                  myBDQR,
                  M_velocityFESpaceETA,
                  M_velocityFESpaceETA,
                  value(density * thickness) * dot (phi_i, phi_j)
        ) >> M_boundaryMass;

        M_boundaryMass->globalAssemble();
        *spcast<MATRIXEPETRA>(convert<BlockMatrix>(mat)->block(0,0)->data()) += *M_boundaryMass;


    }
    shp<aVector>
    FSIAssemblerFE::
    getRightHandSide(const double& time, const shp<aVector>& sol){

        shp<aVector> rhs(NavierStokesAssemblerFE::getRightHandSide(time, sol));


        shp<aVector> retVec(M_BoundaryStiffness->multiplyByVector(sol));

        this->addForcingTerm(retVec);
        retVec->add(rhs);

        return retVec;
        //return rhs;

    }
    shp<aMatrix>
    FSIAssemblerFE::
    getMass(const double& time, const shp<aVector>& sol)//M_mass modificata nel setup una voltasola
    {
        /*shp<BlockMatrix> M_massCopy(new BlockMatrix(this->M_nComponents,this->M_nComponents));
        M_massCopy->deepCopy(NavierStokesAssemblerFE::
                             getMass(time, sol));
        this->addFSIMassMatrix(M_massCopy);
        this->M_bcManager->apply0DirichletMatrix(*M_massCopy,this->getFESpaceBCs(),this->getComponentBCs(), 1.0);

        printlog(YELLOW, "[FSI]Mass matrix assembled ...\n", this->M_data.getVerbose());
        return M_massCopy;*/
        shp<BlockMatrix> retMat(new BlockMatrix(0,0));
        retMat->deepCopy(this->M_mass);
        if (M_stabilization)
        {
            retMat->add(M_stabilization->getMass(convert<BlockVector>(sol),
                                                 convert<BlockVector>(this->getForcingTerm(time))));

            this->M_bcManager->apply0DirichletMatrix(*retMat, this->getFESpaceBCs(),
                                                     this->getComponentBCs(), 1.0);
        }

        return retMat;


    }


    void
    FSIAssemblerFE::
    getBoundaryStiffnessMatrix() const
    {
        using namespace LifeV;
        using namespace ExpressionAssembly;


            shp<MATRIXEPETRA>  M_b(new MATRIXEPETRA(M_velocityFESpace->map()));
            //M_BoundaryStiffness.reset(new MATRIXEPETRA(M_velocityFESpace->map()));
            printlog(YELLOW, "[NS]Assembling boundary stiffness matrix ...\n", this->M_data.getVerbose());
            LifeV::QuadratureBoundary myBDQR(LifeV::buildTetraBDQR(LifeV::quadRuleTria4pt));

            LifeV::MatrixSmall<3, 3> Eye;
            Eye *= 0.0;
            Eye[0][0] = 1;
            Eye[1][1] = 1;
            Eye[2][2] = 1;
            unsigned int wallFlag = M_data("structure/flag", 10);

            double dt = M_data("time_discretization/dt", 0.01);

            integrate(boundary(M_velocityFESpaceETA->mesh(), wallFlag),
                      myBDQR,
                      M_velocityFESpaceETA,
                      M_velocityFESpaceETA,
                      -dt * M_TMAlgorithm->getCoefficientExtrapolation() * (
                              2 * M_lameII *
                              0.5 * dot(
                                      (grad(phi_j) - grad(phi_j) * outerProduct(Nface, Nface))
                                      + transpose(grad(phi_j) - grad(phi_j) * outerProduct(Nface, Nface)),
                                      (grad(phi_i) - grad(phi_i) * outerProduct(Nface, Nface))) +
                              M_lameI *
                              dot(value(Eye), (grad(phi_j) - grad(phi_j) * outerProduct(Nface, Nface))) *
                              dot(value(Eye), (grad(phi_i) - grad(phi_i) * outerProduct(Nface, Nface))))
            ) >> M_b;

            M_b->globalAssemble();

        M_BoundaryStiffness->setBlock(0,0,wrap(M_b));
           // *spcast<MATRIXEPETRA>(M_BoundaryStiffness)->block(0,0)->data()) = *M_b;//ERROR HERE
           this->M_bcManager->apply0DirichletMatrix(*spcast<BlockMatrix>(M_BoundaryStiffness), this->getFESpaceBCs(),
                                                     this->getComponentBCs(), 0.0);



    }
    void
    FSIAssemblerFE::
    addForcingTerm(const shp<aVector> & rhs  ) const
    {

        shp<aVector> rhsDisplacement(M_TMAlgorithm->getPreviousContribution()); //guardare segno


        rhs->add(spcast<BlockMatrix>(M_BoundaryStiffness)->multiplyByVector(rhsDisplacement));

    }

    shp<aMatrix>
    FSIAssemblerFE::
    getJacobianRightHandSide(const double& time,
                             const shp<aVector>& sol)
    {
        shp<aMatrix> retMat = NavierStokesAssemblerFE::getJacobianRightHandSide(time, sol);

        retMat->add(M_BoundaryStiffness);


        this->M_bcManager->apply0DirichletMatrix(*spcast<BlockMatrix>(retMat), this->getFESpaceBCs(),this->getComponentBCs(), 0.0);

        return retMat;
    }

    void
    FSIAssemblerFE::
    postProcess(const double& time, const shp<aVector>& sol)  {
        NavierStokesAssemblerFE::postProcess(time, sol);
        printlog(YELLOW, "[FSI] Updating displacements field ...\n",
                 this->M_data.getVerbose());

        double  dt = this->M_data("time_discretization/dt", 0.01);

        shp<BlockVector> Displacement(new BlockVector(this->M_nComponents));
        M_TMAlgorithm->advanceDisp(dt, convert<BlockVector>(sol));
        printlog(YELLOW, "[FSI] Updating displacements field ...\n",
                 this->M_data.getVerbose());

        Displacement->deepCopy(M_TMAlgorithm->advanceDisp(dt, convert<BlockVector>(sol)));
        printlog(YELLOW, "[FSI] Updating displacements field ...\n",
                 this->M_data.getVerbose());

        M_TMAlgorithm->shiftSolutions(Displacement);
        //*M_displacementExporter = *spcast<VECTOREPETRA>(Displacement->block(0)->data());

    }
    void
    FSIAssemblerFE::
    setExporterDisplacement() {

        printlog(GREEN, "[FSI] Preparing to export solution...\n",
                 this->M_data.getVerbose());
        M_displacementExporter.reset(new VECTOREPETRA(M_velocityFESpace->map(),M_exporter->mapType()));
        M_exporter->addVariable(LifeV::ExporterData<MESH >::VectorField, "displacement", M_velocityFESpace,M_displacementExporter, 0.0);
    }
}
