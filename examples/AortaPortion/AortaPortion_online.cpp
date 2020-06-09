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
#include <redma/solver/problem/DataContainer.hpp>
#include <redma/solver/problem/ProblemRB.hpp>
#include <chrono>
#include <thread>

using namespace RedMA;

int main(int argc, char **argv)
{
    #ifdef HAVE_MPI
    MPI_Init (nullptr, nullptr);
    EPETRACOMM comm (new Epetra_MpiComm(MPI_COMM_WORLD));
    #else
    EPETRACOMM comm(new Epetra_SerialComm());
    #endif

    Epetra_SerialDenseSolver solver;
    DENSEMATRIX matrix(2,2);
    DENSEVECTOR vector(2);
    DENSEVECTOR res(2);

    matrix(0,0) = 2;
    matrix(0,1) = 1;
    matrix(1,0) = -0.5;
    matrix(1,1) = 10;

    DENSEMATRIX matrix2(2,2);
    matrix2(0,0) = 2;
    matrix2(0,1) = 1;
    matrix2(1,0) = -0.5;
    matrix2(1,1) = 10;

    vector(0) = 0;
    vector(1) = 0.5;

    solver.SetMatrix(matrix);
    solver.SetVectors(res, vector);
    solver.Solve();

    DENSEVECTOR residual(2);
    matrix2.Multiply(false, res, residual);
    residual.Scale(-1);
    residual += vector;
    std::cout << "residual = " << residual.Norm2() << std::endl;

    // Chrono chrono;
    // chrono.start();
    //
    // std::string msg = "Starting chrono\n";
    // printlog(MAGENTA, msg, true);
    //
    // DataContainer data;
    // data.setDatafile("datafiles/data");
    // data.setVerbose(comm->MyPID() == 0);
    // data.finalize();
    //
    // ProblemRB rbProblem(data, comm);
    // rbProblem.solve();
    //
    // msg = "Total time =  ";
    // msg += std::to_string(chrono.diff());
    // msg += " seconds\n";
    // printlog(MAGENTA, msg, true);

    return 0;
}
