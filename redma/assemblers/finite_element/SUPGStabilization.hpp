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

#ifndef SUPGSTABILIZATION_HPP
#define SUPGSTABILIZATION_HPP

#include <redma/assemblers/finite_element/NavierStokesStabilization.hpp>

namespace RedMA
{

class SUPGStabilization : public NavierStokesStabilization
{
public:
    SUPGStabilization(const DataContainer& data,
                      SHP(FESPACE) fespaceVelocity,
                      SHP(FESPACE) fespacePressure,
                      SHP(ETFESPACE3) etfespaceVelocity,
                      SHP(ETFESPACE1) etfespacePressure);

    virtual BlockMatrix getMass(const BlockVector& sol,
                                const BlockVector& rhs) override;

    virtual BlockMatrix getMassJac(const BlockVector& sol,
                                   const BlockVector& rhs) override;

    virtual BlockMatrix getJac(const BlockVector& sol,
                               const BlockVector& rhs) override;

    virtual BlockVector getResidual(const BlockVector& sol,
                                    const BlockVector& rhs) override;

};

}  // namespace RedMA

#endif  // SUPGSTABILIZATION_HPP