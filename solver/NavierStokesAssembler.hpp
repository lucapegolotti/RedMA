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

#ifndef NAVIERSTOKESASSEMBLER_HPP
#define NAVIERSTOKESASSEMBLER_HPP

#include <AbstractAssembler.hpp>

#include <lifev/eta/fem/ETFESpace.hpp>

namespace RedMA
{

class NavierStokesAssembler : public AbstractAssembler
{
protected:
    typedef LifeV::ETFESpace<Mesh, Map, 3, 3>   ETFESpaceVelocity;
    typedef std::shared_ptr<ETFESpaceVelocity>  ETFESpaceVelocityPtr;
    typedef LifeV::ETFESpace<Mesh, Map, 3, 1>   ETFESpacePressure;
    typedef std::shared_ptr<ETFESpacePressure>  ETFESpacePressurePtr;


public:
    NavierStokesAssembler(const GetPot& datafile, commPtr_Type comm,
                          const TreeNodePtr& treeNode, bool verbose = false);

    void setup();

protected:
    void assembleConstantMatrices();

    void assembleStiffnessMatrix();
    void assembleDivergenceMatrix();
    void assembleMassMatrix();

    FESpacePtr              M_velocityFESpace;
    FESpacePtr              M_pressureFESpace;
    ETFESpaceVelocityPtr    M_velocityFESpaceETA;
    ETFESpacePressurePtr    M_pressureFESpaceETA;

    MatrixPtr               M_A;
    MatrixPtr               M_B;
    MatrixPtr               M_M;
};

}  // namespace RedMA

#endif  // NAVIERSTOKESASSEMBLER_HPP
