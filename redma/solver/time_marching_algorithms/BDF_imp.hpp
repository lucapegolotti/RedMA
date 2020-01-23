namespace RedMA
{

template <class DataType>
BDF<DataType>::
BDF(const GetPot& datafile) :
  aTimeMarchingAlgorithm<DataType>(datafile),
  M_order(datafile("time_discretization/order",2))
{
    setup();
}

template <class DataType>
void
BDF<DataType>::
setup()
{
    M_coefficients.reserve(M_order);
    M_prevSolutions.reserve(M_order);

    if (M_order == 1)
    {
        M_coefficients[0] = 1;
        M_rhsCoeff;
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

template <class DataType>
DataType
BDF<DataType>::
advance(const double& time, double& dt, SHP(aAssembler) assembler)
{

}

}
