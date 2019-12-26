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

#ifndef ABSTRACTFACTORY_H
#define ABSTRACTFACTORY_H

#include <AbstractAssembler.hpp>
#include <NavierStokesAssembler.hpp>
#include <PseudoFSIAssembler.hpp>

namespace RedMA
{

template <typename... Args>
std::shared_ptr<AbstractAssembler> AssemblerFactory(const GetPot& datafile, Args&&... args)
{
    std::shared_ptr<AbstractAssembler> newAssembler;
    std::string assemblertype = datafile("assembler/type","navierstokes");

    if (!std::strcmp(assemblertype,"navierstokes"))
        newAssembler = new NavierStokesAssembler(args);
    else if (!std::strcmp(assemblertype,"pseusofsi"))
        newAssembler = new PseudoFSIAssembler(args);

    return newAssembler;
}

} // namespace RedMA

#endif
