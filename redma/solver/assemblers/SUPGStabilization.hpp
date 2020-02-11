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

#include <redma/RedMA.hpp>

#include <redma/solver/array/BlockVector.hpp>
#include <redma/solver/array/BlockMatrix.hpp>
#include <redma/solver/array/VectorEp.hpp>
#include <redma/solver/array/MatrixEp.hpp>
#include <redma/utils/Exception.hpp>
#include <redma/solver/problem/DataContainer.hpp>

#include <lifev/eta/utils/Functions.hpp>
#include <lifev/eta/expression/Integrate.hpp>

#include <functional>

namespace RedMA
{

class SUPGStabilization
{
public:
    SUPGStabilization(const DataContainer& data,
                      SHP(FESPACE) fespaceVelocity,
                      SHP(FESPACE) fespacePressure,
                      SHP(ETFESPACE3) etfespaceVelocity,
                      SHP(ETFESPACE1) etfespacePressure);

    void setDensityAndViscosity(const double& density, const double& viscosity);

    BlockMatrix<MatrixEp> getMass(const BlockVector<VectorEp>& sol,
                                  const BlockVector<VectorEp>& rhs);

    BlockMatrix<MatrixEp> getMassJac(const BlockVector<VectorEp>& sol,
                                     const BlockVector<VectorEp>& rhs);

    BlockMatrix<MatrixEp> getJac(const BlockVector<VectorEp>& sol,
                                 const BlockVector<VectorEp>& rhs);

    // BlockMatrix<MatrixEp> assembleMass(const BlockVector<VectorEp>& sol);

    BlockVector<VectorEp> getResidual(const BlockVector<VectorEp>& sol,
                                      const BlockVector<VectorEp>& rhs);

private:
    double                          M_density;
    double                          M_viscosity;
    unsigned int                    M_timeOrder;
    unsigned int                    M_velocityOrder;
    SHP(FESPACE)                    M_velocityFESpace;
    SHP(FESPACE)                    M_pressureFESpace;
    SHP(ETFESPACE3)                 M_velocityFESpaceETA;
    SHP(ETFESPACE1)                 M_pressureFESpaceETA;
    BlockMatrix<MatrixEp>           M_jac;
    BlockMatrix<MatrixEp>           M_massJac;
    BlockMatrix<MatrixEp>           M_mass;
    double                          M_C_I;
    double                          M_dt;
};

}  // namespace RedMA

#endif  // SUPGSTABILIZATION_HPP
