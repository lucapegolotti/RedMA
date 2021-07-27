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
#include <redma/assemblers/finite_element/StokesAssemblerFE.hpp>
#include <redma/problem/DataContainer.hpp>

#include <tests/AtomicTest.hpp>

#include <time.h>

using namespace RedMA;

EPETRACOMM generateComm()
{
    #ifdef HAVE_MPI
    EPETRACOMM comm (new Epetra_MpiComm(MPI_COMM_WORLD));
    #else
    EPETRACOMM comm(new Epetra_SerialComm ());
    #endif
    return comm;
}

shp<StokesAssemblerFE> generateAssembler()
{
    EPETRACOMM comm = generateComm();

    DataContainer data;
    data.setDatafile("datafiles/data");

    shp<Tube> tube(new Tube(comm, "coarse", false, 1, 1, false));
    tube->readMesh("../../../meshes/");
    shp<TreeNode> treeNode(new TreeNode(tube, 0));
    shp<StokesAssemblerFE> assembler(new StokesAssemblerFE(data, treeNode));

    return assembler;
}

int constructor()
{
    shp<StokesAssemblerFE> assembler = generateAssembler();

    return SUCCESS;
}

int setup()
{
    shp<StokesAssemblerFE> assembler = generateAssembler();
    assembler->setup();

    return SUCCESS;
}

int checkNormMatrices()
{
    shp<StokesAssemblerFE> assembler = generateAssembler();
    assembler->setup();

    auto matrices = assembler->getMatrices();
    double normMass = matrices[0]->normInf();
    double normStiffness = matrices[1]->normInf();
    double normDivergence = matrices[2]->normInf();

    int status = SUCCESS;

    status |= (std::abs(normMass - 1.0) > 1e-5);
    status |= (std::abs(normStiffness - 0.699198) > 1e-5);
    status |= (std::abs(normDivergence - 0.506977) > 1e-5);

    return status;
}

int checkNormNorms()
{
    shp<StokesAssemblerFE> assembler = generateAssembler();
    assembler->setup();

    auto normVelocity = assembler->getNorm(0,false)->normInf();
    auto normPressure = assembler->getNorm(1,false)->normInf();

    std::cout << "=====" << std::endl << std::flush;
    std::cout << normVelocity << std::endl << std::flush;
    std::cout << normPressure << std::endl << std::flush;

    int status = SUCCESS;

    status |= (std::abs(normVelocity - 8.82383) > 1e-4);
    status |= (std::abs(normPressure - 0.0164942) > 1e-4);

    return status;
}

int main()
{
    MPI_Init(nullptr, nullptr);
    srand(time(NULL));

    int status = SUCCESS;

    status |= AtomicTest("Constructor", &constructor).run();
    status |= AtomicTest("Setup", &setup).run();
    status |= AtomicTest("CheckNormMatrices", &checkNormMatrices).run();
    status |= AtomicTest("CheckNormNorms", &checkNormNorms).run();

    return status;
}
