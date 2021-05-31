//
// Created by micol on 20/04/2021.
//

#include "FSIAssemblerRB.hpp"


namespace RedMA
{


    FSIAssemblerRB::
    FSIAssemblerRB(const DataContainer& data,
                            shp<TreeNode> treeNode) :
            NavierStokesAssemblerRB(data,treeNode),M_BoundaryStiffnessreduced(new BlockMatrix(2,2))
    ,M_TMAlgorithm(TimeMarchingAlgorithmFactory(data)),M_FSIAssembler(new FSIAssemblerFE(data, treeNode))
    {

        M_name = "FSIAssemblerRB";
        M_exactJacobian = M_data("rb/online/exactjacobian",0);

    }

    shp<aVector>
    FSIAssemblerRB::
    getRightHandSide(const double& time,
                     const shp<aVector>& sol)
    {
        Chrono chrono;
        chrono.start();

        auto solBlck = convert<BlockVector>(sol);

        std::string msg = "[FSIAssemblerRB] computing rhs ...";
        printlog(YELLOW, msg, this->M_data.getVerbose());

        shp<BlockVector> retVec = convert<BlockVector>(NavierStokesAssemblerRB::getRightHandSide(time,sol));

        shp<aVector> rhs(M_BoundaryStiffnessreduced->multiplyByVector(sol));
        double  dt = this->M_data("time_discretization/dt", 0.01);
        rhs->multiplyByScalar(-dt * M_TMAlgorithm->getCoefficientExtrapolation());
        this->addForcingTermreduced(rhs);
        retVec->add(rhs);

        //getBCManager()->apply0DirichletBCs(*spcast<BlockVector>(retVec),getFESpaceBCs(),getComponentBCs());
        msg = "done, in ";
        msg += std::to_string(chrono.diff());
        msg += " seconds\n";
        printlog(YELLOW, msg, this->M_data.getVerbose());

        return retVec;
    }
    void
    FSIAssemblerRB::
    addForcingTermreduced(const shp<aVector> & rhs  ) const
    {

        shp<aVector> rhsDisplacement(M_TMAlgorithm->getPreviousContribution());

        rhs->add(spcast<BlockMatrix>(M_BoundaryStiffnessreduced)->multiplyByVector(rhsDisplacement));

    }

    shp<aMatrix>
    FSIAssemblerRB::getJacobianRightHandSide(const double &time, const shp <aVector> &sol)
    {
        shp<aMatrix> jac = NavierStokesAssemblerRB::getJacobianRightHandSide(time,sol);


            shp<BlockMatrix> retMatcopy(new BlockMatrix(0,0));
            retMatcopy->deepCopy(this->M_BoundaryStiffnessreduced);

            double  dt = this->M_data("time_discretization/dt", 0.01);
            retMatcopy->multiplyByScalar(-dt * M_TMAlgorithm->getCoefficientExtrapolation());
            jac->add(retMatcopy);
            //getBCManager()->apply0DirichletMatrix(*spcast<BlockMatrix>(jac),getFESpaceBCs(),getComponentBCs(),0.0);



        return jac;
    }

    void
    FSIAssemblerRB::
    RBsetup()
    {
        printlog(YELLOW, "[FSIAssemblerRB] assembling and projecting boundary stiffness matrix \t", M_data.getVerbose());
        NavierStokesAssemblerRB::RBsetup();
        M_FSIAssembler->setup();
        M_TMAlgorithm->setComm(M_comm);

        M_TMAlgorithm->setup(this->getZeroVector());

        unsigned int id=M_treeNode->M_ID;
        auto BoundaryStiffness=M_FSIAssembler->returnBoundaryStiffness();
        M_BoundaryStiffnessreduced->setBlock(0,0,M_bases->matrixProject(BoundaryStiffness->block(0,0),0,0,id));
        this->addFSIMassMatrix(M_reducedMass);

        //this->M_bcManager->apply0DirichletMatrix(*M_reducedMass,this->getFESpaceBCs(),this->getComponentBCs(), 1.0);


    }

    void
    FSIAssemblerRB::
    postProcess(const double& time, const shp<aVector>& sol)  {
        NavierStokesAssemblerRB::postProcess(time, sol);
        printlog(YELLOW, "[FSIAssemblerRB] Updating displacements field ...\n",
                 this->M_data.getVerbose());
        unsigned int id=M_treeNode->M_ID;
        double  dt = this->M_data("time_discretization/dt", 0.01);

        shp<BlockVector> Displacement(new BlockVector(this->M_nComponents));



        Displacement->deepCopy(M_TMAlgorithm->advanceDisp(dt, convert<BlockVector>(sol)));
        shp<VECTOREPETRA> exportedDisplacement = M_bases->reconstructFEFunction(Displacement->block(0), 0, id);
        *exportedDisplacement*= *(this->M_FSIAssembler->getBoundaryIndicator());
        M_TMAlgorithm->shiftSolutions(Displacement);
        //*M_displacementExporter = *spcast<VECTOREPETRA>(Displacement->block(0)->data());
        //*M_displacementExporter *= (*M_boundaryIndicator);


    }

    void
    FSIAssemblerRB::
    addFSIMassMatrix(shp<aMatrix> mat)
    {
        Chrono chrono;
        chrono.start();

        std::string msg = "[FSIAssemblerRB] computing boundary mass ...";
        printlog(YELLOW, msg, this->M_data.getVerbose());
        using namespace LifeV;
        using namespace ExpressionAssembly;


        LifeV::QuadratureBoundary myBDQR(LifeV::buildTetraBDQR(LifeV::quadRuleTria4pt)); //not sure if updated, there is another way
        shp<MATRIXEPETRA> M_boundaryMass(new MATRIXEPETRA(M_feStokesAssembler->getFEspace(0)->map()));
        shp<ETFESPACE3> velocityFESpaceETA = M_feStokesAssembler->getVelocityETFEspace();
        double density = M_feStokesAssembler->getDensity();
        double thickness = M_data("structure/thickness", 0.1);
        unsigned int wallFlag = M_data("structure/flag", 10);
        integrate(boundary(velocityFESpaceETA->mesh(), wallFlag),
                  myBDQR,
                  velocityFESpaceETA,
                  velocityFESpaceETA,
                  value(density * thickness) * dot (phi_i, phi_j)
        ) >> M_boundaryMass;

        M_boundaryMass->globalAssemble();

        shp<SparseMatrix> matrixWrap(new SparseMatrix());
        matrixWrap->setMatrix(M_boundaryMass);

        unsigned int id = M_treeNode->M_ID;
        auto projectedMat = M_bases->matrixProject(matrixWrap, 0, 0, id);

        convert<BlockMatrix>(mat)->block(0,0)->add(projectedMat);

        msg = "done, in ";
        msg += std::to_string(chrono.diff());
        msg += " seconds\n";
        printlog(YELLOW, msg, this->M_data.getVerbose());



    }

}
