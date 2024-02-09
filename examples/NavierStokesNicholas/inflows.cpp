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

double inflow_bypass(const double t, const std::vector<double> params, const double Tramp)
{
    BSpline spline;
    spline.knots = {0.000000,  0.000000,  0.000000,  0.000000,  0.013009,  0.019696,
                    0.026383,  0.033678,  0.040365,  0.043404,  0.048875,  0.053131,
                    0.059210,  0.064073,  0.070152,  0.075623,  0.079271,  0.083526,
                    0.093252,  0.102979,  0.107842,  0.112705,  0.120608,  0.130334,
                    0.142492,  0.152827,  0.163161,  0.168632,  0.176535,  0.183830,
                    0.191125,  0.197204,  0.201459,  0.206930,  0.214833,  0.220304,
                    0.224559,  0.230638,  0.235502,  0.240365,  0.244620,  0.247660,
                    0.251307,  0.255562,  0.260426,  0.266505,  0.270152,  0.275623,
                    0.279271,  0.285957,  0.292644,  0.299331,  0.304802,  0.310274,
                    0.315137,  0.320000,  0.323040,  0.324255,  0.326687,  0.327903,
                    0.330942,  0.332766,  0.333982,  0.335198,  0.337021,  0.340669,
                    0.341884,  0.344924,  0.349179,  0.357082,  0.364377,  0.374711,
                    0.383222,  0.389301,  0.397812,  0.403283,  0.407538,  0.412401,
                    0.414225,  0.416049,  0.419088,  0.422736,  0.424559,  0.431246,
                    0.440973,  0.443404,  0.446444,  0.453739,  0.465897,  0.481702,
                    0.495076,  0.502371,  0.513921,  0.522432,  0.533982,  0.547356,
                    0.563769,  0.580182,  0.588693,  0.597204,  0.607538,  0.619088,
                    0.633070,  0.644012,  0.651307,  0.661641,  0.673799,  0.685957,
                    0.698723,  0.716960,  0.733982,  0.743100,  0.753435,  0.764985,
                    0.800000,  0.800000,  0.800000,  0.800000};
    spline.controlPoints = {0.000000,  0.107424,  0.270066,  0.686590,  1.079219,  1.529738,
                            1.665662,  2.479341,  2.526228,  3.203902,  3.330118,  3.823002,
                            4.134922,  4.419553,  4.866427,  5.221333,  5.486763,  5.489934,
                            5.871342,  6.121987,  6.415745,  6.376537,  6.730727,  6.422683,
                            6.752654,  6.440733,  6.332223,  6.216715,  5.968089,  5.894305,
                            5.600864,  5.428963,  5.369481,  5.061707,  4.794320,  4.624618,
                            4.313969,  4.128478,  3.923578,  3.642798,  3.488203,  3.259535,
                            2.999609,  2.896841,  2.565957,  2.385971,  1.959699,  1.810823,
                            1.575021,  1.394963,  1.243907,  1.005325,  0.896644,  0.603720,
                            0.632064,  0.022120,  0.143219, -0.469150, -0.431882, -0.740658,
                            -1.010876, -1.553619, -1.760727, -1.723078, -2.329966, -2.367022,
                            -2.828083, -2.655143, -2.557138, -2.541135, -2.543657, -2.546730,
                            -2.530823, -2.369128, -2.139345, -2.180916, -1.869986, -1.511433,
                            -1.222019, -1.160536, -0.412524, -0.584778, -0.679712, -0.082808,
                            0.411558,  0.472179,  0.526822,  0.530207,  0.607694,  0.662179,
                            0.784511,  0.610707,  0.861281,  0.720515,  0.747265,  0.519615,
                            0.806465,  0.930635,  0.743332,  0.720410,  0.765768,  0.603564,
                            0.547469,  0.661133,  0.604920,  0.617585,  0.628546,  0.546654,
                            0.432488,  0.558382,  0.369634,  0.407125,  0.174742,  0.000000};
    std::transform(spline.controlPoints.begin(), spline.controlPoints.end(), spline.controlPoints.begin(),
                   [](auto& c){return -1.0 * c;});


    spline.degree = 3;

    return evaluateBSpline(spline, t);

}

double splineBasisFunction(unsigned int i, unsigned int p, const std::vector<double>& knots, double t) {
    if (p == 0)
        return (knots[i] <= t && t < knots[i + 1]) ? 1.0 : 0.0;

    double denom1 = knots[i + p] - knots[i];
    double denom2 = knots[i + p + 1] - knots[i + 1];

    double alpha1 = (denom1 == 0) ? 0.0 : (t - knots[i]) / denom1;
    double alpha2 = (denom2 == 0) ? 0.0 : (knots[i + p + 1] - t) / denom2;

    double retVal = (alpha1 * splineBasisFunction(i, p - 1, knots, t) +
                     alpha2 * splineBasisFunction(i + 1, p - 1, knots, t));

    return retVal;
}


double evaluateBSpline(const BSpline& spline, double t) {
    double result = 0.0;

    for (unsigned int i = 0; i < spline.knots.size() - spline.degree - 1; ++i)
        result += splineBasisFunction(i, spline.degree, spline.knots, t) * spline.controlPoints[i];

    return result;
}
