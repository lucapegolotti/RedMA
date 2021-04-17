#include "BDF.hpp"

namespace RedMA {

    BDF::
    BDF(const DataContainer &data) :
            aTimeMarchingAlgorithm(data),
            M_order(data("time_discretization/order", 2)) {
        M_useExtrapolation = this->M_data("time_discretization/use_extrapolation", 0);
        // if we set this we save the evaluation of the residual at the end of resolution
        // (important for rb method)
        if (M_useExtrapolation)
            this->M_systemSolver.isLinearProblem();
    }

    BDF::
    BDF(const DataContainer &data, shp <FunProvider> funProvider) :
            aTimeMarchingAlgorithm(data, funProvider),
            M_order(data("time_discretization/order", 2)) {
        setup(this->M_funProvider->getZeroVector());

        M_useExtrapolation = this->M_data("time_discretization/use_extrapolation", 0);

        if (M_useExtrapolation)
            this->M_systemSolver.isLinearProblem();
    }


    void
    BDF::
    setup(const shp <aVector> &zeroVector) {
        M_coefficients.reserve(M_order);

        for (unsigned int i = 0; i < M_order; i++) {
            shp <BlockVector> zeroVectorCopy(new BlockVector(0));
            zeroVectorCopy->deepCopy(zeroVector);
            M_prevSolutions.push_back(zeroVectorCopy);
        }

        if (M_order == 1) {
            M_coefficients[0] = -1.0;
            M_rhsCoeff = 1.0;
        } else if (M_order == 2) {
            M_coefficients[0] = -4.0 / 3.0;
            M_coefficients[1] = 1.0 / 3.0;
            M_rhsCoeff = 2.0 / 3.0;
        } else if (M_order == 3) {
            M_coefficients[0] = -18.0 / 11.0;
            M_coefficients[1] = 9.0 / 11.0;
            M_coefficients[2] = -2.0 / 11.0;
            M_rhsCoeff = 6.0 / 11.0;
        } else {
            throw new Exception("BDF scheme of requested order not implemented");
        }
    }

    shp <aVector>
    BDF::
    computeExtrapolatedSolution() {
        shp <BlockVector> extrapolatedSolution(new BlockVector(0));
        extrapolatedSolution->deepCopy(M_prevSolutions[0]);
        // std::cout << "extrapolatedSolution open?" << extrapolatedSolution->isOpen() << std::endl << std::flush;
        // std::cout << "extrapolatedSolution open?" << extrapolatedSolution->isOpen() << std::endl << std::flush;
        if (M_order == 1) {

        } else if (M_order == 2) {
            // std::cout << "extrapolatedSolution open?" << extrapolatedSolution->isOpen() << std::endl << std::flush;
            extrapolatedSolution->multiplyByScalar(2);
            M_prevSolutions[1]->multiplyByScalar(-1);
            extrapolatedSolution->add(M_prevSolutions[1]);
            M_prevSolutions[1]->multiplyByScalar(-1);
        } else if (M_order == 3) {
            extrapolatedSolution->multiplyByScalar(3);
            // M_prevSolutions[1]->close();
            M_prevSolutions[1]->multiplyByScalar(-3);
            extrapolatedSolution->add(M_prevSolutions[1]);
            M_prevSolutions[1]->multiplyByScalar(-1.0 / 3.0);
            // M_prevSolutions[2]->close();
            extrapolatedSolution->add(M_prevSolutions[2]);
        } else
            throw new Exception("BDF scheme of requested order not implemented");

        return extrapolatedSolution;
    }

    shp <aVector>
    BDF::
    advance(const double &time, double &dt, int &status) {
        typedef shp <aVector> BV;
        typedef shp <aMatrix> BM;

        // we set the initial guess equal to the last solution
        // keep in mind that this MUST be a hard copy
        BV initialGuess = computeExtrapolatedSolution();

        this->M_funProvider->applyDirichletBCs(time + dt, initialGuess);

        FunctionFunctor <BV, BV> fct(
                [this, time, dt](BV sol) {
                    BM mass(this->M_funProvider->getMass(time + dt, sol));
                    if (M_useExtrapolation)
                        this->M_funProvider->setExtrapolatedSolution(computeExtrapolatedSolution());

                    BV f(this->M_funProvider->getRightHandSide(time + dt, sol));

                    BV prevContribution(new BlockVector(0));
                    unsigned int count = 0;
                    for (BV vec : M_prevSolutions) {
                        BV vecCopy(new BlockVector(0));
                        vecCopy->deepCopy(M_prevSolutions[count]);
                        vecCopy->multiplyByScalar(M_coefficients[count]);
                        prevContribution->add(vecCopy);
                        count++;
                    }
                    prevContribution->add(sol);

                    BV retVec(new BlockVector(0));
                    retVec->deepCopy(mass->multiplyByVector(prevContribution));

                    f->multiplyByScalar(-1. * M_rhsCoeff * dt);
                    retVec->add(f);
                    // the previous solution satisfies the boundary conditions so we search
                    // for an increment with 0bcs
                    this->M_funProvider->apply0DirichletBCs(retVec);
                    return retVec;
                });

        FunctionFunctor <BV, BM> jac(
                [this, time, dt](BV sol) {
                    // here the choice of hard copies is compulsory
                    BM retMat(new BlockMatrix(1, 1));

                    if (M_useExtrapolation)
                        this->M_funProvider->setExtrapolatedSolution(computeExtrapolatedSolution());
                    retMat->deepCopy(this->M_funProvider->getJacobianRightHandSide(time + dt, sol));
                    retMat->multiplyByScalar(-1. * M_rhsCoeff * dt);
                    retMat->add(this->M_funProvider->getMass(time + dt, sol));
                    retMat->add(this->M_funProvider->getMassJacobian(time + dt, sol));

                    // this is the part relative to the previous steps
                    unsigned int count = 0;
                    for (BV vec : M_prevSolutions) {
                        auto massjac = this->M_funProvider->getMassJacobian(time + dt, vec);
                        massjac->multiplyByScalar(M_coefficients[count]);
                        retMat->add(massjac);
                        count++;
                    }
                    return retMat;
                });
        BV sol = this->M_systemSolver.solve(fct, jac, initialGuess, status);

        if (status != 0)
            throw new Exception("Solver has not converged!");

        std::vector <SolverStatistics> statistics = this->M_systemSolver.getSolverStatistics();
        this->dumpSolverStatistics(statistics, time + dt);

        return sol;
    }

    shp <aVector>
    BDF::
    computeDerivative(const shp <aVector> &solnp1, double &dt) {
        typedef shp <BlockVector> BV;

        BV retVec(new BlockVector(0));
        retVec->deepCopy(solnp1);

        unsigned int count = 0;
        for (BV vec : M_prevSolutions) {
            BV vecCopy(new BlockVector(0));
            vecCopy->deepCopy(vec);
            vecCopy->multiplyByScalar(M_coefficients[count]);
            retVec->add(vecCopy);
            count++;
        }

        retVec->multiplyByScalar(1.0 / (dt * M_rhsCoeff));
        return retVec;
    }

    void
    BDF::
    shiftSolutions(const shp <aVector> &sol) {
        // shift solutions
        std::vector <shp<BlockVector>> newPrevSolutions(M_order);
        newPrevSolutions[0].reset(static_cast<BlockVector *>(sol->clone()));

        for (unsigned int i = 0; i < M_order - 1; i++)
            newPrevSolutions[i + 1].reset(new BlockVector(*M_prevSolutions[i]));

        M_prevSolutions = newPrevSolutions;
    }

    double
    BDF::
    getCoefficientExtrapolation() {
        return M_rhsCoeff;
    }


    shp <aVector>
    BDF::
    getPreviousContribution() {
        typedef shp <aVector> BV;
        BV prevContribution(new BlockVector(0));
        unsigned int count = 0;
        for (BV vec : M_prevSolutions) {
            BV vecCopy(new BlockVector(0));
            vecCopy->deepCopy(M_prevSolutions[count]);
            vecCopy->multiplyByScalar(M_coefficients[count]);
            prevContribution->add(vecCopy);
            count++;
        }
        return prevContribution;
    }
    shp<aVector>
    BDF::
    advanceDisp(const double &dt, const shp<BlockVector> &sol) {
        typedef shp <BlockVector> BV;
        shp<BlockVector> retVec(new BlockVector(2));


        retVec->deepCopy(sol);
        retVec->multiplyByScalar(dt * M_rhsCoeff);
        shp<BlockVector> oldDisp(new BlockVector(2));

        unsigned int count = 0;
        for (BV vec : M_prevSolutions) {
            BV vecCopy(new BlockVector(2));//distrvect

            (vecCopy)->deepCopy(M_prevSolutions[count]);
            vecCopy->multiplyByScalar(-M_coefficients[count]);
            oldDisp->add(vecCopy);
            count++;
        }

        retVec->block(0)->add(oldDisp->block(0));
        return retVec;
    }

    shp<aVector>
    BDF::getPreviousSolution(){
        return   M_prevSolutions.end()[-2];
    }



}
