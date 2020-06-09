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
#include <time.h>       /* time */

using namespace RedMA;

int main(int argc, char **argv)
{
    srand(1234);

    #ifdef HAVE_MPI
    MPI_Init (nullptr, nullptr);
    EPETRACOMM comm (new Epetra_MpiComm(MPI_COMM_WORLD));
    #else
    EPETRACOMM comm(new Epetra_SerialComm());
    #endif

    unsigned int N = atoi(argv[1]);

    Epetra_SerialDenseSolver solver;
    DENSEMATRIX matrix(N,N);
    DENSEMATRIX matrix2(N,N);
    DENSEVECTOR vector(N);
    DENSEVECTOR res(N);

    for (unsigned int i = 0; i < N; i++)
    {
        int num1 = rand() % 100 - 50;
        vector(i) = num1;
        for (unsigned int j = 0; j < N; j++)
        {
            int num2 = rand() % 100 - 50;
            matrix(i,j) = num2;
            matrix2(i,j) = num2;
        }
    }

    solver.SolveToRefinedSolution(true);
    solver.SetMatrix(matrix);
    solver.SetVectors(res, vector);
    solver.Solve();

    DENSEVECTOR residual(N);
    matrix2.Multiply(false, res, residual);
    residual.Scale(-1);
    residual += vector;
    std::cout << "factored = " << solver.Factored() << std::endl;
    std::cout << "factored matrix = " << solver.FactoredMatrix()->NormInf() << std::endl;
    std::cout << "residual = " << residual.Norm2() << std::endl;
    std::cout << "ferr = " << *solver.FERR() << std::endl;
    std::cout << "berr = " << *solver.BERR() << std::endl;

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
