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

#ifndef STOKESASSEMBLER_HPP
#define STOKESASSEMBLER_HPP

#include <redma/solver/assemblers/aAssembler.hpp>

namespace RedMA
{

class StokesAssembler : public aAssembler
{
public:
    StokesAssembler(const GetPot& datafile);

    virtual void exportSolution(const double& t) override;

    virtual void postProcess() override;

protected:
};

}

#endif // STOKESASSEMBLER_HPP
