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

#ifndef COMMUNICATORDISTRIBUTOR_HPP
#define COMMUNICATORDISTRIBUTOR_HPP

#include <redma/geometry/GeometryParser.hpp>

namespace RedMA
{

class CommunicatorsDistributor
{
    typedef std::shared_ptr<Epetra_Comm>  commPtr_Type;

public:
    CommunicatorsDistributor(const GetPot& datafile, commPtr_Type comm);

    void loadBalancing(std::map<unsigned int, unsigned int> npointsmeshes);

    void createCommunicators();

    inline std::vector<commPtr_Type> getComms() const {return M_communicators;}

    inline std::map<unsigned int, unsigned int> getProcessMap() const {return M_processMap;}

private:
    GetPot                                                  M_datafile;
    commPtr_Type                                            M_comm;
    std::shared_ptr<GeometryParser>                         M_parser;
    std::vector<commPtr_Type>                               M_communicators;
    TreeStructure                                           M_tree;
    // key = building block id, value = process id
    std::map<unsigned int, unsigned int>                    M_processMap;
};

}  // namespace RedMA

#endif  // COMMUNICATORDISTRIBUTOR_HPP
