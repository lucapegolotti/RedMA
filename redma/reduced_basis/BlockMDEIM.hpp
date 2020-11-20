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

#ifndef BLOCKMDEIM_HPP
#define BLOCKMDEIM_HPP

#include <redma/RedMA.hpp>
#include <redma/problem/DataContainer.hpp>
#include <redma/reduced_basis/RBBases.hpp>
#include <redma/array/BlockMatrix.hpp>
#include <redma/array/BlockVector.hpp>

// #include <rb/reduced_basis/rbSolver/ProperOrthogonalDecomposition.hpp>
// #include <rb/reduced_basis/util/EpetraArrayUtils.hpp>
#include <redma/reduced_basis/MDEIM.hpp>

#include <filesystem>

namespace RedMA
{

class BlockMDEIM
{
    typedef boost::numeric::ublas::matrix<shp<MDEIM>> GridMDEIM;

public:
    BlockMDEIM();

    void setDataContainer(const DataContainer& dataContainer);

    void setComm(EPETRACOMM comm);

    void addSnapshot(BlockMatrix newSnapshot);

    void performMDEIM(std::string outdir);

    void initialize(BlockMatrix reducedMatrix);

    void checkOnline(BlockMatrix completeMat, BlockMatrix reducedMat);

    void checkOnSnapshots();

    void dumpMDEIMs(std::string dir, std::string prefix = "blockmdeim");

    void projectMDEIMs();

    void loadMDEIM(std::string dir);

    void setRBBases(shp<RBBases> bases);

    inline void setMatrixIndex(const unsigned int& index) {M_matIndex = index;}

    void resize(unsigned int rows, unsigned int cols);

    BlockMatrix assembleMatrix(BlockMatrix reducedMatrix);

    BlockMatrix assembleProjectedMatrix(BlockMatrix reducedMatrix);

    void prepareOnline();

    void setFESpace(shp<FESPACE> fespace, const unsigned int& index);

    void setRangeMap(shp<MAPEPETRA> map, const unsigned int& index);

    void setDomainMap(shp<MAPEPETRA> map, const unsigned int& index);

    inline shp<FESPACE> getFESpace(const unsigned int& index) {return M_fespaces[index];}

    inline shp<MAPEPETRA> getRangeMap(const unsigned int& index) {return M_rangeMaps[index];}

    inline shp<MAPEPETRA> getDomainMap(const unsigned int& index) {return M_domainMaps[index];}

    inline unsigned int getNumRows() const {return M_nRows;}

    inline unsigned int getNumCols() const {return M_nCols;}

    inline unsigned int getIndex() const {return M_matIndex;}

    inline BlockMDEIMStructure& getMDEIMStructure() {return M_structures;};

private:

    GridMDEIM                                   M_mdeims;
    unsigned int                                M_matIndex;
    EPETRACOMM                                  M_comm;
    DataContainer                               M_data;
    unsigned int                                M_nRows;
    unsigned int                                M_nCols;
    BlockMDEIMStructure                         M_structures;
    shp<RBBases>                                M_bases;
    std::vector<shp<FESPACE>>                   M_fespaces;
    std::vector<shp<MAPEPETRA>>                 M_rangeMaps;
    std::vector<shp<MAPEPETRA>>                 M_domainMaps;
};

}  // namespace RedMA

#endif  // BLOCKMDEIM_HPP
