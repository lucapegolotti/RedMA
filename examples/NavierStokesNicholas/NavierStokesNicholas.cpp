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

#include <redma/RedMA.hpp>
#include <redma/problem/GlobalProblem.hpp>
#include <redma/problem/DataContainer.hpp>

#include <cmath>
#include<Eigen/Dense>
#include <thread>

#include <redma/reduced_basis/SnapshotsSampler.hpp>

using namespace RedMA;

double inflow(const double t, const std::vector<double> params, const double T)
{
    return (1-cos(2*M_PI*t/T) + params[1]*sin(2*M_PI*params[0]*t/T));
}


double inflow_systolic(const double t, const std::vector<double> params)
{
    double Tramp = 0.05;

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


double inflow_heartbeat(const double t, const std::vector<double> params)
{
    double Tramp = 0.05;

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


int main(int argc, char **argv)
{

    std::mt19937_64 eng{std::random_device{}()};
    std::uniform_int_distribution<> dist{1, 20};
    std::this_thread::sleep_for(std::chrono::seconds{dist(eng)});

    Chrono chrono;
    chrono.start();

    std::string msg = "Starting chrono \n";
    printlog(MAGENTA, msg, true);
    
    #ifdef HAVE_MPI
    MPI_Init (nullptr, nullptr);
    EPETRACOMM comm (new Epetra_MpiComm(MPI_COMM_WORLD));
    #else
    EPETRACOMM comm(new Epetra_SerialComm());
    #endif

    printlog(MAGENTA,"Starting snapshots generation", true);
    DataContainer data;
    data.setDatafile("datafiles/data_fem");
    data.setVerbose(comm->MyPID() == 0);

    unsigned int Nstart = 0;
    if (argc > 1)
        Nstart = std::atoi(argv[1]);

    double T = data("time_discretization/T", 1.0);

    if (std::strcmp(data("rb/offline/snapshots/param_type", "inflow").c_str(), "inflow"))
        throw new Exception("This test case handles only 'inflow' parametrization!");

    SnapshotsSampler sampler(data, comm);
    if (!std::strcmp(data("rb/offline/snapshots/inflow_type", "inflow").c_str(), "default"))
        sampler.setInflow([T](const double t, const std::vector<double> params){return inflow(t, params, T);});
    else if (!std::strcmp(data("rb/offline/snapshots/inflow_type", "inflow").c_str(), "systolic"))
        sampler.setInflow(inflow_systolic);
    else if (!std::strcmp(data("rb/offline/snapshots/inflow_type", "inflow").c_str(), "heartbeat"))
        sampler.setInflow(inflow_heartbeat);
    else
        throw new Exception("Unrecognized type of inflow parametrization! "
                            "Available types: {default, systolic, heartbeat}.");

    sampler.takeSnapshots(Nstart);

    // To solve the RB problem in RedMA  --> set datafile accordingly!
    /*data.setInletBC(inflow2, 0);
    data.finalize();
    GlobalProblem rbProblem(data, comm);
    rbProblem.solve();*/

    // To compute the convective term, given a solution
    /*GlobalProblem femProblem(data, comm);
    femProblem.validateRBConvectiveTerm();*/

    msg = "Total time =  ";
    msg += std::to_string(chrono.diff());
    msg += " seconds\n";
    printlog(MAGENTA, msg, true);

    return 0;
}
