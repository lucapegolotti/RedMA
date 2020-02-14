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

#include <redma/solver/assemblers/NavierStokesStabilization.hpp>


#define ALPHA 1e-5

namespace RedMA
{

class HFStabilization : public NavierStokesStabilization
{
public:
    HFStabilization(const DataContainer& data,
                    SHP(FESPACE) fespaceVelocity,
                    SHP(FESPACE) fespacePressure,
                    SHP(ETFESPACE3) etfespaceVelocity,
                    SHP(ETFESPACE1) etfespacePressure);

    virtual BlockMatrix<MatrixEp> getMass(const BlockVector<VectorEp>& sol,
                                          const BlockVector<VectorEp>& rhs) override;

    virtual BlockMatrix<MatrixEp> getMassJac(const BlockVector<VectorEp>& sol,
                                             const BlockVector<VectorEp>& rhs) override;

    virtual BlockMatrix<MatrixEp> getJac(const BlockVector<VectorEp>& sol,
                                         const BlockVector<VectorEp>& rhs) override;

    virtual BlockVector<VectorEp> getResidual(const BlockVector<VectorEp>& sol,
                                              const BlockVector<VectorEp>& rhs) override;

};

}  // namespace RedMA

#endif  // SUPGSTABILIZATION_HPP
