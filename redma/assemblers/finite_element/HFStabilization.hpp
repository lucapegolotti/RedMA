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

#ifndef HFSTABILIZATION_HPP
#define HFSTABILIZATION_HPP

#include <redma/assemblers/finite_element/NavierStokesStabilization.hpp>


#define ALPHA 1e-5

namespace RedMA
{

class HFStabilization : public NavierStokesStabilization
{
public:
    HFStabilization(const DataContainer& data,
                    shp<FESPACE> fespaceVelocity,
                    shp<FESPACE> fespacePressure,
                    shp<ETFESPACE3> etfespaceVelocity,
                    shp<ETFESPACE1> etfespacePressure);

    virtual shp<BlockMatrix> getMass(shp<BlockVector> sol,
                                     shp<BlockVector> rhs) override;

    virtual shp<BlockMatrix> getMassJac(shp<BlockVector> sol,
                                        shp<BlockVector> rhs) override;

    virtual shp<BlockMatrix> getJac(shp<BlockVector> sol,
                                    shp<BlockVector> rhs) override;

    virtual shp<BlockVector> getResidual(shp<BlockVector> sol,
                                         shp<BlockVector> rhs) override;

};

}  // namespace RedMA

#endif  // SUPGSTABILIZATION_HPP
