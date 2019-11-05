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
#include "WomersleySolution.hpp"

#include <lifev/core/filter/GetPot.hpp>

#include "ud_functions.hpp"

using namespace RedMA;

double maxLaw(double t)
{
    return LifeV::pressureWomerTime(t,0,0,0,0);
}

double maxLawDt(double t)
{
    return LifeV::pressuredtWomerTime(t,0,0,0,0);
}

int main(int argc, char **argv)
{
    #ifdef HAVE_MPI
    MPI_Init (nullptr, nullptr);
    std::shared_ptr<Epetra_Comm> comm (new Epetra_MpiComm(MPI_COMM_WORLD));
    #else
    std::shared_ptr<Epetra_Comm> comm(new Epetra_SerialComm());
    #endif

    GetPot datafile("datafiles/data");
    LifeV::Womersley::setParamsFromGetPot(datafile);
    bool verbose = comm->MyPID() == 0;

    std::vector<std::string> refinements;
    // refinements.push_back("h0.80");
    // refinements.push_back("h0.70");
    // refinements.push_back("h0.60");
    // refinements.push_back("h0.50");
    // refinements.push_back("h0.40");
    // refinements.push_back("h0.30");
    // refinements.push_back("h0.25");
    // refinements.push_back("h0.20");
    refinements.push_back("h0.10");

    for (auto it = refinements.begin(); it != refinements.end(); it++)
    {
        std::string nameTree("datafiles/tree_");
        nameTree += *it + std::string(".xml");

        datafile.set("geometric_structure/xmlfile",
                     nameTree.c_str());

        AbstractFunctor* womerlseySolution = new WomersleySolution;
        GlobalSolver<NavierStokesAssembler> gs(datafile, comm, verbose,
                                               womerlseySolution);
        gs.setExportErrors("errors" + *it + ".txt");
        gs.printMeshSize("meshSizes" + *it + ".txt");

        gs.setLawInflow(std::function<double(double)>(maxLaw));
        gs.setLawDtInflow(std::function<double(double)>(maxLawDt));
        gs.solve();

        delete womerlseySolution;
    }

    return 0;
}
