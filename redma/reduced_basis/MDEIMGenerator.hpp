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

#ifndef MDEIMGENERATOR_HPP
#define MDEIMGENERATOR_HPP

#include <redma/RedMA.hpp>
#include <redma/solver/problem/DataContainer.hpp>
#include <redma/solver/problem/ProblemFEM.hpp>
#include <redma/geometry/GeometryPrinter.hpp>

#include <redma/reduced_basis/BlockMDEIM.hpp>

#include <boost/filesystem.hpp>

namespace RedMA
{

class MDEIMGenerator
{
public:

    MDEIMGenerator(const DataContainer& data, EPETRACOMM comm);

    void generateMDEIM();

private:
    void takeMatricesSnapshots();

    void performMDEIM();

    void checkMDEIM();

    void projectMDEIM();

    void dumpMDEIMstructures();

    std::string getMeshName(std::string mapkey);

    DataContainer                                     M_data;
    EPETRACOMM                                        M_comm;
    // one MDEIM instance for every matrix of the problem (e.g. stiffness, mass, ecc)
    std::map<std::string, std::vector<BlockMDEIM>>    M_primalBlockMDEIMsMap;
    // first mdeim: inlet, then outlets in order
    // std::map<std::string, std::vector<BlockMDEIM>>    M_dualBlockTMDEIMsMap;
    // std::map<std::string, std::vector<BlockMDEIM>>    M_dualBlockMDEIMsMap;

};

}  // namespace RedMA

#endif  // MDEIMGENERATOR_HPP
