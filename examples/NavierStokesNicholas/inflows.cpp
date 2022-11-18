#include "inflows.hpp"

using namespace RedMA;

double inflow(const double t, const std::vector<double> params, const double T)
{
    return (1-cos(2*M_PI*t/T) + params[1]*sin(2*M_PI*params[0]*t/T));
}

double inflow_periodic(const double t, const std::vector<double> params, const double T, const double Tramp)
{
    if (Tramp < 1e-8)
        throw new Exception("Ramp time must be strictly positive!");

    if (t<0)
        return 5.0 * (1 - cos((t+Tramp) * M_PI / Tramp));
    else
        return 10.0 + std::abs(params[1]*sin(2*M_PI*params[0]*t/T));
}


double inflow_systolic(const double t, const std::vector<double> params, const double Tramp)
{
    if (Tramp < 1e-8)
        throw new Exception("Ramp time must be strictly positive!");

    // reference values, computed from measured inflow
    double V0_ref = 1.541;
    double TM_ref = 0.13375;
    double VM_ref = 14.161;
    double Ts_ref = 0.3075;
    double Tm_ref = 0.375;
    double Vm_ref = 0.626;

    double V0 = V0_ref * (1.0 + params[1]);    // initial flow
    double TM = TM_ref * (1.0 + params[0]);    // time of systolic peak
    double VM = VM_ref * (1.0 + params[2]);    // peak systolic flow
    double Ts = Ts_ref * (1.0 + params[0]);    // systolic time
    double Tm = Tm_ref * (1.0 + params[0]);    // time to min flow
    double Vm = Vm_ref * (1.0 + params[3]);    // min flow

    Eigen::Matrix<double, 8, 8> matrix_sys;
    matrix_sys << 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0,
    pow(TM,7), pow(TM,6), pow(TM,5), pow(TM,4), pow(TM,3), pow(TM,2), TM, 1.0,
    pow(Ts,7), pow(Ts,6), pow(Ts,5), pow(Ts,4), pow(Ts,3), pow(Ts,2), Ts, 1.0,
    pow(Tm,7), pow(Tm,6), pow(Tm,5), pow(Tm,4), pow(Tm,3), pow(Tm,2), Tm, 1.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0,
    7.*pow(TM,6), 6.*pow(TM,5), 5.*pow(TM,4), 4.*pow(TM,3), 3.*pow(TM,2), 2.*TM, 1.0, 0.0,
    7.*pow(Tm,6), 6.*pow(Tm,5), 5.*pow(Tm,4), 4.*pow(Tm,3), 3.*pow(Tm,2), 2.*Tm, 1.0, 0.0,
    42.*pow(Tm,5), 30.*pow(Tm,4), 20.*pow(Tm,3), 12.*pow(Tm,2), 6.*Tm, 2.0, 0.0, 0.0;

    Eigen::Matrix<double, 8, 1> vector_sys;
    vector_sys << V0, VM, V0, Vm, 0.0, 0.0, 0.0, 0.0;

    Eigen::Matrix<double, 8, 1> a_sys = matrix_sys.colPivHouseholderQr().solve(vector_sys);
    FunctionFunctor<double, double> systolic_flow(
            [a_sys](double t)
            {
                return a_sys[0]*pow(t,7) + a_sys[1]*pow(t,6) + a_sys[2]*pow(t,5) + a_sys[3]*pow(t,4) +
                a_sys[4]*pow(t,3) + a_sys[5]*pow(t,2) + a_sys[6]*pow(t,1) + a_sys[7]*pow(t,0);
            });

    if (t<0)
        return (V0/2) * (1 - cos((t+Tramp) * M_PI / Tramp));
    else
        return systolic_flow(t);
}


double inflow_heartbeat(const double t, const std::vector<double> params, const double Tramp)
{
    if (Tramp < 1e-8)
        throw new Exception("Ramp time must be strictly positive!");

    // reference values, computed from measured inflow
    double V0_ref = 1.541;
    double TM_ref = 0.13375;
    double VM_ref = 14.161;
    double Ts_ref = 0.3075;
    double Tm_ref = 0.375;
    double Vm_ref = 0.626;
    double TMd_ref = 0.63375;
    double VMd_ref = 2.092;
    double Tf_ref = 0.75;
    double Vf_ref = 1.527;

    double V0 = V0_ref * (1.0 + params[1]);    // initial flow
    double TM = TM_ref * (1.0 + params[0]);    // time of systolic peak
    double VM = VM_ref * (1.0 + params[2]);    // peak systolic flow
    double Ts = Ts_ref * (1.0 + params[0]);    // systolic time
    double Tm = Tm_ref * (1.0 + params[0]);    // time to min flow
    double Vm = Vm_ref * (1.0 + params[3]);    // min flow
    double TMd = TMd_ref;                      // time of diastolic peak
    double VMd = VMd_ref * (1.0 + params[4]);  // peak diastolic flow
    double Tf = Tf_ref;                        // final time
    double Vf = Vf_ref * (1.0 + params[1]);    // flow at final time

    Eigen::Matrix<double, 8, 8> matrix_sys;
    matrix_sys << 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0,
    pow(TM,7), pow(TM,6), pow(TM,5), pow(TM,4), pow(TM,3), pow(TM,2), TM, 1.0,
    pow(Ts,7), pow(Ts,6), pow(Ts,5), pow(Ts,4), pow(Ts,3), pow(Ts,2), Ts, 1.0,
    pow(Tm,7), pow(Tm,6), pow(Tm,5), pow(Tm,4), pow(Tm,3), pow(Tm,2), Tm, 1.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0,
    7.*pow(TM,6), 6.*pow(TM,5), 5.*pow(TM,4), 4.*pow(TM,3), 3.*pow(TM,2), 2.*TM, 1.0, 0.0,
    7.*pow(Tm,6), 6.*pow(Tm,5), 5.*pow(Tm,4), 4.*pow(Tm,3), 3.*pow(Tm,2), 2.*Tm, 1.0, 0.0,
    42.*pow(Tm,5), 30.*pow(Tm,4), 20.*pow(Tm,3), 12.*pow(Tm,2), 6.*Tm, 2.0, 0.0, 0.0;

    Eigen::Matrix<double, 8, 1> vector_sys;
    vector_sys << V0, VM, V0, Vm, 0.0, 0.0, 0.0, 0.0;

    Eigen::Matrix<double, 8, 1> a_sys = matrix_sys.colPivHouseholderQr().solve(vector_sys);
    FunctionFunctor<double, double> systolic_flow(
            [a_sys](double t)
            {
                return a_sys[0]*pow(t,7) + a_sys[1]*pow(t,6) + a_sys[2]*pow(t,5) + a_sys[3]*pow(t,4) +
                a_sys[4]*pow(t,3) + a_sys[5]*pow(t,2) + a_sys[6]*pow(t,1) + a_sys[7]*pow(t,0);
            });


    double Td = Tm_ref;
    if (Tm < Tm_ref)
        Td = Tm;

    double Vd = systolic_flow(Td);
    double Vdp = (systolic_flow(Td)- systolic_flow(Td-0.001)) / 0.001;

    Eigen::Matrix<double, 5, 5> matrix_dia;
    matrix_dia << pow(Td,4), pow(Td,3), pow(Td,2), Td, 1.0,
    4.*pow(Td,3), 3.*pow(Td,2), 2.*Td, 1.0, 0.0,
    pow(TMd,4), pow(TMd,3), pow(TMd,2), TMd, 1.0,
    4.*pow(TMd,3), 3.*pow(TMd,2), 2.*TMd, 1.0, 0.0,
    pow(Tf,4), pow(Tf,3), pow(Tf,2), Tf, 1.0;

    Eigen::Matrix<double, 5, 1> vector_dia;
    vector_dia << Vd, Vdp, VMd, 0.0, Vf;

    Eigen::Matrix<double, 5, 1> a_dia = matrix_dia.colPivHouseholderQr().solve(vector_dia);
    FunctionFunctor<double, double> diastolic_flow(
            [a_dia](double t)
            {
                return a_dia[0]*pow(t,4) + a_dia[1]*pow(t,3) + a_dia[2]*pow(t,2) +
                a_dia[3]*pow(t,1) + a_dia[4]*pow(t,0);
            });

    if (t<0)
        return (V0/2) * (1 - cos((t+Tramp) * M_PI / Tramp));
    else if (t<Td)
        return systolic_flow(t);
    else
        return diastolic_flow(t);
}