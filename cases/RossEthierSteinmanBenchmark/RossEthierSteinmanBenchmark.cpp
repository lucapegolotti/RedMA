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

#include <Epetra_ConfigDefs.h>
#ifdef EPETRA_MPI
#include <mpi.h>
#include <Epetra_MpiComm.h>
#else
#include <Epetra_SerialComm.h>
#endif

#include <functional>

#include <GlobalSolver.hpp>
#include <NavierStokesAssembler.hpp>
#include <PseudoFSIAssembler.hpp>
#include "RossEthierSteinmanSolution.hpp"
#include <fstream>
#include <lifev/core/filter/GetPot.hpp>

using namespace RedMA;

int main(int argc, char **argv)
{
    #ifdef HAVE_MPI
    MPI_Init (nullptr, nullptr);
    std::shared_ptr<Epetra_Comm> comm (new Epetra_MpiComm(MPI_COMM_WORLD));
    #else
    std::shared_ptr<Epetra_Comm> comm(new Epetra_SerialComm());
    #endif

    GetPot datafile("datafiles/data");
    LifeV::RossEthierSteinmanDec::setParamsFromGetPot(datafile);
    bool verbose = comm->MyPID() == 0;

    std::vector<std::string> basisFunctions;
    basisFunctions.push_back("chebyshev");
    basisFunctions.push_back("zernike");

    std::vector<std::string> refinements;
    refinements.push_back("h0.80");
    refinements.push_back("h0.70");
    refinements.push_back("h0.60");
    refinements.push_back("h0.50");
    refinements.push_back("h0.40");
    refinements.push_back("h0.30");
    refinements.push_back("h0.25");
    refinements.push_back("h0.20");
    refinements.push_back("h0.15");
    refinements.push_back("h0.13");

    for (int nMax = 0; nMax < 3; nMax++)
    {
        std::string nMaxString(std::to_string(nMax));
        datafile.set("coupling/nMax", nMaxString.c_str());

        for (auto itBfs = basisFunctions.begin(); itBfs != basisFunctions.end(); itBfs++)
        {
            datafile.set("coupling/type", itBfs->c_str());

            for (auto it = refinements.begin(); it != refinements.end(); it++)
            {
                std::string nameTree("datafiles/tree_");
                nameTree += *it + std::string(".xml");

                datafile.set("geometric_structure/xmlfile", nameTree.c_str());

                AbstractFunctor* RESSolution = new RossEthierSteinmanSolution;
                GlobalSolver<NavierStokesAssembler> gs(datafile, comm, verbose,
                                                       RESSolution);

                gs.setForcingFunction(LifeV::RossEthierSteinmanDec::f,
                                      LifeV::RossEthierSteinmanDec::f_dt);

                std::string errorFile("errors");
                errorFile += "_";
                errorFile += *itBfs;
                errorFile += "_";
                errorFile += nMaxString;
                errorFile += "_";
                errorFile += *it;
                errorFile += ".txt";
                std::ifstream f(errorFile.c_str());
                if (!f.good())
                {
                    gs.setExportErrors(errorFile);
                    gs.printMeshSize("meshSizes" + *it + ".txt");
                    gs.solve();
                }
                delete RESSolution;
            }
        }
    }

    return 0;
}
