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

    bool verbose = comm->MyPID() == 0;

    if (argc < 4)
    {
        std::cout << "Insufficient number of arguments!" << std::endl;
        exit(1);
    }

    std::string bf = argv[1];
    std::string nMax = argv[2];
    std::string ref = "h" + std::string(argv[3]);

    GetPot datafile("datafiles/data");
    LifeV::RossEthierSteinmanDec::setParamsFromGetPot(datafile);
    datafile.set("coupling/nMax", nMax.c_str());
    datafile.set("coupling/type", bf.c_str());

    std::string nameTree("datafiles/tree_");
    nameTree += ref + std::string(".xml");

    datafile.set("geometric_structure/xmlfile", nameTree.c_str());

    AbstractFunctor* RESSolution = new RossEthierSteinmanSolution;
    GlobalSolver<NavierStokesAssembler> gs(datafile, comm, verbose,
                                           RESSolution);

    gs.setForcingFunction(LifeV::RossEthierSteinmanDec::f,
                          LifeV::RossEthierSteinmanDec::f_dt);

    std::string errorFile("errors");
    errorFile += "_";
    errorFile += bf;
    errorFile += "_";
    errorFile += nMax;
    errorFile += "_";
    errorFile += ref;
    errorFile += ".txt";
    std::ifstream f(errorFile.c_str());

    gs.setExportErrors(errorFile);
    gs.printMeshSize("meshSizes" + ref + ".txt");
    gs.solve();
    delete RESSolution;

    return 0;
}
