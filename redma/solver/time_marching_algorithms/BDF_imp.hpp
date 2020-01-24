namespace RedMA
{

template <class InVectorType, class InMatrixType>
BDF<InVectorType, InMatrixType>::
BDF(const GetPot& datafile) :
  aTimeMarchingAlgorithm<InVectorType, InMatrixType>(datafile),
  M_order(datafile("time_discretization/order",2))
{
    setup();
}

template <class InVectorType, class InMatrixType>
void
BDF<InVectorType, InMatrixType>::
setup()
{
    M_coefficients.reserve(M_order);
    M_prevSolutions.reserve(M_order);

    if (M_order == 1)
    {
        M_coefficients[0] = 1.0;
        M_rhsCoeff = 1.0;
    }
    else if (M_order == 2)
    {
        M_coefficients[0] = -4.0/3.0;
        M_coefficients[1] = 1.0/3.0;
        M_rhsCoeff = 2.0/3.0;
    }
    else if (M_order == 3)
    {
        M_coefficients[0] = -18.0/11.0;
        M_coefficients[1] = 9.0/11.0;
        M_coefficients[2] = -2.0/11.0;
        M_rhsCoeff = 6.0/11.0;
    }
    else
    {
        throw new Exception("BDF scheme of requested order not implemented");
    }
}

template <class InVectorType, class InMatrixType>
BlockVector<InVectorType>
BDF<InVectorType, InMatrixType>::
advance(const double& time, double& dt,
        SHP((aAssembler<InVectorType, InMatrixType>)) assembler)
{
    typedef BlockVector<InVectorType>               BV;
    typedef BlockMatrix<InMatrixType>               BM;

    FunctionFunctor<BV,BV> fct(
        [this,time,dt,assembler](BV sol)
    {
        BM mass = assembler->getMass(time+dt, sol);
        BV f = assembler->getRightHandSide(time+dt, sol);
        BV prevContribution;

        unsigned int count = 0;
        for (BV vec : M_prevSolutions)
        {
            prevContribution += vec * M_coefficients[count];
            count++;
        }

        BV retVec;
        retVec = mass * (sol + prevContribution);
        f *= (-1. * M_rhsCoeff * dt);
        retVec += f;

        return retVec;
    });

    FunctionFunctor<BV,BM> jac(
        [this,time,dt,assembler](BV sol)
    {
        BM retMat = assembler->getMass(time+dt, sol);
        BM jacRhs;

        jacRhs = assembler->getJacobianRightHandSide(time+dt, sol);
        jacRhs *= (-1. * M_rhsCoeff * dt);

        retMat += jacRhs;
        return retMat;
    });

    return M_systemSolver.solve(fct,jac);
}

}
