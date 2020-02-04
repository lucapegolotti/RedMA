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
#include <redma/solver/problem/ProblemFEM.hpp>
#include <redma/solver/problem/DataContainer.hpp>
#include <redma/geometry/CommunicatorsDistributor.hpp>

using namespace RedMA;

double inflow(double t)
{
    return 10;
}

double inflowDt(double t)
{
    return 0;
}

void solveProblem(EPETRACOMM comm, std::string innerprec, std::string approxschur,
                  double intol, int numblocks, int nmax)
{
    DataContainer data;
    data.setDatafile("datafiles/data");
    data.setInflow(inflow);
    data.setInflowDt(inflowDt);
    data.setVerbose(comm->MyPID() == 0);

    data.setValue("geometric_structure/maxnumblocks", numblocks);
    data.setValue("preconditioner/inner", innerprec);
    data.setValue("preconditioner/approxschur", approxschur);
    data.setValue("preconditioner/innertol", intol);
    data.setValue("coupling/nMax", nmax);
    std::string outfile = "statistics";
    outfile += "mnb" + std::to_string(numblocks);
    outfile += "in" + innerprec;
    outfile += "as" + approxschur;
    outfile += "intol" + std::to_string(intol);
    outfile += "nMax" + std::to_string(nmax);
    outfile += ".txt";
    std::cout << "Writing " << outfile << std::endl;
    data.setValue("time_discretization/solverstatistics", outfile);

    ProblemFEM femProblem(data, comm);

    femProblem.solve();
}

int main(int argc, char **argv)
{
    #ifdef HAVE_MPI
    MPI_Init (nullptr, nullptr);
    EPETRACOMM comm (new Epetra_MpiComm(MPI_COMM_WORLD));
    #else
    EPETRACOMM comm(new Epetra_SerialComm());
    #endif

    std::vector<std::string> innerPrec(2); innerPrec[0] = "SIMPLE"; innerPrec[1] = "exact";
    std::vector<std::string> approxSchur(2); approxSchur[0] = "SIMPLE"; approxSchur[1] = "exact";
    std::vector<double> intol(3); intol[0] = 0.5; intol[1] = 0.1; intol[2] = 0.01;

    DataContainer data;
    data.setDatafile("datafiles/data");

    CommunicatorsDistributor commD(data.getDatafile(), comm);

    // for (auto inp : innerPrec)
    // {
    //     for (auto aps : approxSchur)
    //     {
    //         if (!std::strcmp(inp.c_str(),"exact") || !std::strcmp(aps.c_str(),"exact"))
    //         {
    //             for (auto tol : intol)
    //             {
    //                 for (int numblocks = 1; numblocks < 40; )
    //                 {
    //                     for (int nmax = 0; nmax < 3; nmax++)
    //                         solveProblem(comm,inp,aps,tol,numblocks,nmax);
    //                     numblocks = numblocks + 3;
    //                 }
    //             }
    //         }
    //         else
    //         {
    //             for (int numblocks = 1; numblocks < 40; )
    //             {
    //                 for (int nmax = 0; nmax < 3; nmax++)
    //                     solveProblem(comm,inp,aps,0.5,numblocks,nmax);
    //                 numblocks = numblocks + 3;
    //             }
    //         }
    //     }
    // }

    #ifdef HAVE_MPI
    MPI_Finalize();
    #endif
    return ( EXIT_SUCCESS );

    return 0;
}
