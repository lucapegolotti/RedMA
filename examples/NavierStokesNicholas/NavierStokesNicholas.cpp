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

double inflow(double t, std::vector<double> params)
{
    double T = 0.3;
    double K = 1;

    return K * (1-cos(2*M_PI*t/T) + params[1]*sin(2*M_PI*params[0]*t/T));
}

double inflow_systolic(double t, std::vector<double> params)
{
    // reference values, computed from measured inflow
    double V0_ref = 1.541;
    double TM_ref = 0.13375;
    double VM_ref = 14.161;
    double Ts_ref = 0.3075;
    double Tm_ref = 0.375;
    double Vm_ref = 0.626;

    double V0 = V0_ref * (1.0 + params[1]);  // initial flow
    double TM = TM_ref * (1.0 + params[0]);  // time of systolic peak
    double VM = VM_ref * (1.0 + params[2]);  // peak systolic flow
    double Ts = Ts_ref * (1.0 + params[0]);  // systolic time
    double Tm = Tm_ref * (1.0 + params[0]);  // time to min flow
    double Vm = Vm_ref * (1.0 + params[3]);  // min flow

    Eigen::Matrix<double, 8, 8> matrix;
    matrix << 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0,
              pow(TM,7), pow(TM,6), pow(TM,5), pow(TM,4), pow(TM,3), pow(TM,2), TM, 1.0,
              pow(Ts,7), pow(Ts,6), pow(Ts,5), pow(Ts,4), pow(Ts,3), pow(Ts,2), Ts, 1.0,
              pow(Tm,7), pow(Tm,6), pow(Tm,5), pow(Tm,4), pow(Tm,3), pow(Tm,2), Tm, 1.0,
              0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0,
              7.*pow(TM,6), 6.*pow(TM,5), 5.*pow(TM,4), 4.*pow(TM,3), 3.*pow(TM,2), 2.*TM, 1.0, 0.0,
              7.*pow(Tm,6), 6.*pow(Tm,5), 5.*pow(Tm,4), 4.*pow(Tm,3), 3.*pow(Tm,2), 2.*Tm, 1.0, 0.0,
              42.*pow(Tm,5), 30.*pow(Tm,4), 20.*pow(Tm,3), 12.*pow(Tm,2), 6.*Tm, 2.0, 0.0, 0.0;

    Eigen::Matrix<double, 8, 1> vector;
    vector << V0, VM, V0, Vm, 0.0, 0.0, 0.0, 0.0;

    Eigen::Matrix<double, 8, 1> a = matrix.colPivHouseholderQr().solve(vector);

    double Tramp = 0.05;

    if (t<0)
        return (V0/2) * (1 - cos((t+Tramp) * M_PI / Tramp));
    else
        return a[0]*pow(t,7) + a[1]*pow(t,6) + a[2]*pow(t,5) + a[3]*pow(t,4) +
               a[4]*pow(t,3) + a[5]*pow(t,2) + a[6]*pow(t,1) + a[7]*pow(t,0);
}

int main(int argc, char **argv)
{

    std::mt19937_64 eng{std::random_device{}()};
    std::uniform_int_distribution<> dist{1, 20};
    std::this_thread::sleep_for(std::chrono::seconds{dist(eng)});

    Chrono chrono;
    chrono.start();

    std::string msg = "Starting chrono... \n";
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

    SnapshotsSampler sampler(data, comm);
    if (!std::strcmp(data("rb/offline/snapshots/param_type", "inflow").c_str(), "inflow"))
        sampler.setInflow(inflow);
    else if (!std::strcmp(data("rb/offline/snapshots/param_type", "inflow").c_str(), "inflow_systolic"))
        sampler.setInflow(inflow_systolic);
    else
        throw new Exception("Unrecognized type of inflow parametrization! "
                            "Available types: {inflow, inflow_systolic}.");

    sampler.takeSnapshots(Nstart);

    msg = "Total time =  ";
    msg += std::to_string(chrono.diff());
    msg += " seconds\n";
    printlog(MAGENTA, msg, true);

    return 0;
}
