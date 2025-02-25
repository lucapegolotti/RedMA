#include "inflows.hpp"

using namespace RedMA;

double inflow(const double t, const std::vector<double> params, const double T,
              const double scale)
{
    return scale * (1-cos(2*M_PI*t/T) + params[1]*sin(2*M_PI*params[0]*t/T));
}

double inflow_periodic(const double t, const std::vector<double> params, const double T, const double Tramp,
                       const double scale)
{
    if (Tramp < 1e-8)
        throw new Exception("Ramp time must be strictly positive!");

    if (t<0)
        return scale * 5.0 * (1 - cos((t+Tramp) * M_PI / Tramp));
    else
        return scale * 10.0 + std::abs(params[1]*sin(2*M_PI*params[0]*fmod(t, T)/T));
}


double inflow_systolic(const double t, const std::vector<double> params, const double T, const double Tramp,
                       const double scale)
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

    assert (T == Tm_ref);

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
        return scale * (V0/2) * (1 - cos((t+Tramp) * M_PI / Tramp));
    else
        return scale * systolic_flow(fmod(t, T));
}


double inflow_heartbeat(const double t, const std::vector<double> params, const double T, const double Tramp,
                        const double scale)
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

    assert (T == Tf_ref);

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
        return scale * (V0/2) * (1 - cos((t+Tramp) * M_PI / Tramp));
    else
    {
        double tMod = fmod(t, T);
        if (tMod < Td)
            return scale * systolic_flow(fmod(t, T));
        else
            return scale * diastolic_flow(fmod(t, T));
    }
}

double inflow_bypass(const double t, const std::vector<double> params, const double T, const double scale)
{
    assert (T == 0.80);

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
                   [scale](auto& c){return scale * c;});

    spline.degree = 3;

    return evaluateBSpline(spline, fmod(t, T));

}

double outflow_bypass(const double t, const std::vector<double> params, const double T, const double scale)
{
    assert (T == 0.80);

    BSpline spline;
    spline.knots = {0.        , 0.        , 0.        , 0.        , 0.00222222,
                    0.00444444, 0.00666667, 0.00888889, 0.01111111, 0.01333333,
                    0.01555556, 0.01777778, 0.02      , 0.03048327, 0.03866171,
                    0.04460967, 0.04907063, 0.05576208, 0.06171004, 0.07137546,
                    0.07732342, 0.08550186, 0.0929368 , 0.10111524, 0.1070632 ,
                    0.11375465, 0.12342007, 0.13457249, 0.14498141, 0.15167286,
                    0.15390335, 0.15687732, 0.15836431, 0.1598513 , 0.1605948 ,
                    0.16133829, 0.16356877, 0.16431227, 0.16802974, 0.17100372,
                    0.1732342 , 0.1739777 , 0.17546468, 0.17620818, 0.17769517,
                    0.17843866, 0.18066914, 0.18215613, 0.18289963, 0.18513011,
                    0.18959108, 0.19330855, 0.20817844, 0.21412639, 0.22007435,
                    0.22379182, 0.23048327, 0.23717472, 0.2394052 , 0.24312268,
                    0.24535316, 0.24832714, 0.25204461, 0.25576208, 0.26245353,
                    0.2669145 , 0.27286245, 0.27509294, 0.2802974 , 0.28475836,
                    0.28698885, 0.29591078, 0.31078067, 0.32490706, 0.34052045,
                    0.35687732, 0.36802974, 0.37695167, 0.38513011, 0.39405204,
                    0.40520446, 0.41412639, 0.42230483, 0.42750929, 0.44386617,
                    0.45427509, 0.46394052, 0.47063197, 0.47881041, 0.49144981,
                    0.49665428, 0.5063197 , 0.51598513, 0.52565056, 0.53457249,
                    0.54200743, 0.54795539, 0.55687732, 0.56133829, 0.5732342 ,
                    0.58364312, 0.59107807, 0.60371747, 0.6133829 , 0.62453532,
                    0.63568773, 0.64684015, 0.65576208, 0.66245353, 0.67509294,
                    0.68475836, 0.69070632, 0.6936803 , 0.69442379, 0.69591078,
                    0.70185874, 0.7063197 , 0.71003717, 0.71747212, 0.72416357,
                    0.72936803, 0.73085502, 0.73605948, 0.73977695, 0.74126394,
                    0.7464684 , 0.75910781, 0.76802974, 0.77546468, 0.78      ,
                    0.78222222, 0.78444444, 0.78666667, 0.78888889, 0.79111111,
                    0.79333333, 0.79555556, 0.79777778, 0.8       , 0.8       ,
                    0.8       , 0.8 };
    spline.controlPoints = {0.        ,  0.0114031 ,  0.03420929,  0.21153903,  0.47112086,
                            0.7923105 ,  1.13319945,  1.45467187,  1.7131226 ,  1.89469396,
                            1.9358361 ,  2.66445842,  3.00006167,  3.41034204,  4.03832239,
                            4.27416276,  4.83065416,  4.78006153,  5.3534451 ,  5.47558548,
                            5.76167234,  5.72553231,  6.0468728 ,  6.22133741,  5.86138032,
                            6.28198243,  6.06018535,  6.07167588,  6.87074263,  6.89855466,
                            7.46259226,  7.61483051,  8.1247303 ,  8.78240768,  8.63081528,
                            9.40117935,  9.26463161,  8.99007553,  9.11903575,  8.28653217,
                            8.37461547,  7.52550251,  7.62907148,  6.67511634,  6.9008124 ,
                            6.53362503,  5.81516343,  5.81198389,  5.46865751,  4.90407916,
                            5.42525376,  4.79951252,  4.59973213,  3.97788861,  3.60759734,
                            3.65470833,  2.79043679,  2.87292451,  2.20343042,  2.00655357,
                            1.59605413,  1.21879554,  1.10980417,  0.40644433,  0.51974372,
                            -0.59203009, -0.7057029 , -0.90252017, -1.92253192, -1.59956304,
                            -1.79772976, -1.62592415, -1.55993208, -1.43202631, -1.71555109,
                            -1.41927449, -1.20820251, -0.95761819, -1.07167944, -0.81594861,
                            -0.71835637,  0.00508391, -0.32814105, -0.00461953,  0.15338827,
                            0.5319553 ,  0.80974215,  0.69682703,  1.04553653,  1.12543914,
                            1.14329098,  1.0640457 ,  0.81610174,  1.12405914,  0.80334519,
                            0.86028034,  1.35298399,  1.01541468,  0.98609567,  1.32340832,
                            1.36757105,  1.32441407,  1.4668423 ,  1.35005258,  1.33145863,
                            1.28879385,  1.06799196,  1.17147405,  1.33184902,  0.90971089,
                            1.10924362,  0.44565021,  0.15958389,  0.07643125,  0.34714656,
                            0.58679833,  0.73667865,  0.42956946,  0.52633164, -0.14420232,
                            -0.12714337, -0.0688376 , -0.94910189, -0.88223117, -0.9226313 ,
                            -0.64090488, -0.52991735, -0.30170202, -0.3161628 , -0.28245275,
                            -0.24071168, -0.18729726, -0.13100867, -0.07788697, -0.03497508,
                            -0.0056556 , -0.0018852 ,  0. };
    std::transform(spline.controlPoints.begin(), spline.controlPoints.end(), spline.controlPoints.begin(),
                   [scale](auto& c){return scale * c;});

    spline.degree = 3;

    return evaluateBSpline(spline, fmod(t, T));

}


double outpres_bypass(double t, const double T, const double scale)
{
    assert (T == 0.80);

    BSpline spline;
    spline.knots = {0.   , 0.   , 0.   , 0.   , 0.014, 0.022, 0.03 , 0.035, 0.043,
                    0.052, 0.06 , 0.068, 0.08 , 0.095, 0.112, 0.129, 0.142, 0.151,
                    0.162, 0.165, 0.17 , 0.182, 0.192, 0.201, 0.212, 0.222, 0.233,
                    0.24 , 0.246, 0.254, 0.263, 0.267, 0.273, 0.279, 0.285, 0.29 ,
                    0.297, 0.303, 0.309, 0.314, 0.32 , 0.327, 0.334, 0.341, 0.35 ,
                    0.359, 0.367, 0.376, 0.386, 0.393, 0.403, 0.412, 0.421, 0.437,
                    0.451, 0.459, 0.473, 0.492, 0.507, 0.52 , 0.536, 0.552, 0.563,
                    0.577, 0.585, 0.596, 0.611, 0.624, 0.638, 0.648, 0.665, 0.68 ,
                    0.693, 0.698, 0.709, 0.722, 0.735, 0.743, 0.756, 0.764, 0.773,
                    0.783, 0.8  , 0.8  , 0.8  , 0.8  };
    spline.controlPoints =  {    0.        ,   -521.48251419,  -1340.95503648,  -4616.50032124,
                             -5828.74893772,  -9660.39357054, -11484.39558431, -13730.75627326,
                             -16570.45064737, -19537.72183305, -21409.16861739, -22956.53550642,
                             -24464.16230353, -23665.3832951 , -23244.67765976, -21401.20380516,
                             -21477.45045257, -25392.52345192, -27375.05453267, -27904.57118393,
                             -30049.69259646, -28434.15133711, -27762.84212459, -24968.07252322,
                             -24321.4425788 , -21725.50176581, -19649.16548502, -16488.15467324,
                             -16479.54528549, -10548.79362969,  -8099.97573919,  -3658.64276386,
                             -2021.39164656,   4235.45853374,   6753.44664826,  12561.0511864 ,
                             15310.9087777 ,  19439.32596497,  24207.08418647,  26722.26442256,
                             29617.93641242,  30508.07176185,  29371.87691459,  27026.93928186,
                             24290.2745283 ,  21904.58433314,  20586.87582401,  16064.78711981,
                             15081.32256866,  10784.49796136,   7396.62131739,   7437.44207104,
                             5595.46941939,   8339.31179722,   7749.37617829,   7849.34419816,
                             6444.42287558,   4749.34900765,   3378.33382674,   3845.41494956,
                             1298.54164124,   2110.25997049,   -242.98006492,  -1362.49383005,
                             -2413.02333098,  -4796.9233633 ,  -6134.39509991,  -9308.39382937,
                             -10908.18714958, -12333.52115868, -15161.70591118, -12045.95571246,
                             -15257.64112962, -12735.90963429, -14062.24764372,  -9362.80338765,
                             -8887.7111589 ,  -5464.6601181 ,  -2633.5842729 ,  -1975.45284671,
                             -763.24314532,      0.  };

    std::transform(spline.controlPoints.begin(), spline.controlPoints.end(), spline.controlPoints.begin(),
                   [scale](auto& c){return scale * c;});

    spline.degree = 3;

    return evaluateBSpline(spline, fmod(t, T));

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
