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
#include <redma/array/BlockMaps.hpp>

#include <tests/AtomicTest.hpp>

#include <time.h>

using namespace RedMA;


shp<Epetra_Comm> generateComm()
{
    #ifdef HAVE_MPI
    std::shared_ptr<Epetra_Comm> comm (new Epetra_MpiComm(MPI_COMM_WORLD));
    #else
    std::shared_ptr<Epetra_Comm> comm(new Epetra_SerialComm ());
    #endif
    return comm;
}

shp<LifeV::MapEpetra> generateMap(unsigned int size, shp<Epetra_Comm> comm)
{
    shp<LifeV::MapEpetra> map;
    map.reset(new LifeV::MapEpetra(size, size, 0, comm));

    return map;
}

shp<DENSEMATRIX> randmatdense(unsigned int N, unsigned int M)
{
    shp<DENSEMATRIX> mat(new DENSEMATRIX(N,M));

    for (unsigned int i = 0; i < N; i++)
    {
        for (unsigned int j = 0; j < M; j++)
        {
            mat->operator()(i,j) = rand() % 10;
        }
    }
    return mat;
}
shp<MATRIXEPETRA> fromDense(shp<DENSEMATRIX> otherMatrix)
{
    int N = otherMatrix->M();
    int M = otherMatrix->N();

    shp<Epetra_Comm> comm = generateComm();
    shp<LifeV::MapEpetra> rowMap = generateMap(N, comm);
    shp<LifeV::MapEpetra> colMap = generateMap(M, comm);

    shp<MATRIXEPETRA> mat(new MATRIXEPETRA(*rowMap,M,false));

    for (unsigned int i = 0; i < N; i++)
    {
        for (unsigned int j = 0; j < M; j++)
            mat->addToCoefficient(i, j, otherMatrix->operator()(i,j));
    }
    mat->globalAssemble(colMap, rowMap);
    return mat;
}

shp<MATRIXEPETRA> randmatsparse(unsigned int N, unsigned int M)
{
    return fromDense(randmatdense(N,M));
}

shp<BlockMatrix> generateBlockMatrix(unsigned int N, unsigned int M,
                                     std::vector<int> sizesrows,
                                     std::vector<int> sizescols)
{
    shp<BlockMatrix> retMat(new BlockMatrix(N,M));

    for (unsigned int i = 0; i < N; i++)
    {
        for (unsigned int j = 0; j < M; j++)
        {
            retMat->setBlock(i,j,wrap(randmatsparse(sizesrows[i],sizescols[j])));
        }
    }

    return retMat;
}

void fillrowscols(int N, int M, std::vector<int>& rows, std::vector<int>& cols)
{
    for (unsigned int i = 0; i < N; i++)
        rows.push_back(i+1);

    for (unsigned int i = M; i > 0; i--)
        cols.push_back(i);
}

int constructor1()
{
    BlockMaps bmaps;
    return SUCCESS;
}

int constructor2()
{
    unsigned int N = 5;

    std::vector<int> rows;
    std::vector<int> cols;
    for (unsigned int i = 0; i < N; i++)
    {
        rows.push_back(2);
        cols.push_back(2);
    }
    fillrowscols(N,N,rows,cols);

    shp<BlockMatrix> bmat1 = generateBlockMatrix(N, N, rows, cols);

    BlockMaps bmaps(bmat1);

    if (bmaps.M_collapsedMatrix->nRows() != N)
        return FAILURE;

    if (bmaps.M_collapsedMatrix->nCols() != N)
        return FAILURE;

    return SUCCESS;
}

int constructor3()
{
    unsigned int N = 5;

    std::vector<int> rows;
    std::vector<int> cols;
    for (unsigned int i = 0; i < N; i++)
    {
        rows.push_back(2);
        cols.push_back(2);
    }
    fillrowscols(N,N,rows,cols);

    shp<BlockMatrix> bmat1 = generateBlockMatrix(N, N, rows, cols);
    shp<BlockMatrix> bmat2 = generateBlockMatrix(N, N, rows, cols);
    shp<BlockMatrix> bmat3 = generateBlockMatrix(N, N, rows, cols);
    shp<BlockMatrix> bmat4 = generateBlockMatrix(N, N, rows, cols);
    shp<BlockMatrix> bmat(new BlockMatrix(2,2));

    bmat->setBlock(0,0,bmat1);
    bmat->setBlock(0,1,bmat2);
    bmat->setBlock(1,0,bmat3);
    bmat->setBlock(1,1,bmat4);

    BlockMaps bmaps(bmat);

    if (bmaps.M_collapsedMatrix->nRows() != N * 2)
        return FAILURE;

    if (bmaps.M_collapsedMatrix->nCols() != N * 2)
        return FAILURE;

    if (bmat->level() != 2)
        return FAILURE;

    if (bmaps.M_collapsedMatrix->level() != 1)
        return FAILURE;

    return SUCCESS;
}


int main()
{
    MPI_Init(nullptr, nullptr);
    srand(time(NULL));

    int status = SUCCESS;

    status |= AtomicTest("constructor1", &constructor1).run();
    status |= AtomicTest("Constructor2", &constructor2).run();
    status |= AtomicTest("Constructor3", &constructor3).run();

    return status;
}
