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

#ifndef LINEAROPERATOR_HPP
#define LINEAROPERATOR_HPP

#include <redma/RedMA.hpp>
#include <redma/array/BlockVector.hpp>
#include <redma/array/BlockMatrix.hpp>
#include <redma/array/BlockMaps.hpp>

#include <lifev/core/filter/GetPot.hpp>
#include <lifev/core/linear_algebra/LinearOperatorAlgebra.hpp>
#include <lifev/core/linear_algebra/BlockEpetra_MultiVector.hpp>
#include <lifev/core/linear_algebra/BlockEpetra_Map.hpp>
#include <lifev/core/linear_algebra/InvertibleOperator.hpp>

namespace RedMA
{

/*! \brief Linear operator.
 *
 * Some methods in this class have an empty implementation, due to
 * the inheritance from LifeV::Operators::LinearOperatorAlgebra.
 */
class LinearOperator : public LifeV::Operators::LinearOperatorAlgebra
{
    typedef LifeV::Operators::LinearOperatorAlgebra         super;
    typedef shp<aVector>                                    BV;
    typedef shp<aMatrix>                                    BM;

public:
    /*! \brief Constructor.
     *
     * \param matrix Shared pointer to the matrix of the operator.
     * \param maps Shared pointer to the block maps.
     */
    LinearOperator(const BM& matrix,
                   shp<BlockMaps> maps);

    /// Method not implemented.
    virtual int SetUseTranspose(bool UseTranspose) override {}

    /*! \brief Apply the linear operator to a vector.
     *
     * \param X Vector to which the operator must be applied.
     * \param Y Result.
     * \return Return code; 0 if successful.
     */
    virtual int Apply(const super::vector_Type& X,
                      super::vector_Type& Y) const override;

    /// Method not implemented.
    virtual int ApplyInverse(const super::vector_Type& X,
                             super::vector_Type& Y) const override {return -1;}

    /// Method not implemented.
    virtual double NormInf() const override {return -1;}

    /// Method not implemented.
    virtual const char * Label() const override {}

    /// Method not implemented.
    virtual bool UseTranspose() const override {return false;}

    /// Method not implemented.
    virtual bool HasNormInf() const override {return false;}

    /*! \brief Getter for the MPI Communicator.
     *
     * \return Shared pointer to the MPI Communicator.
     */
    virtual const super::comm_Type& Comm() const override {return *M_comm;};

    /// Method not implemented.
    virtual const super::map_Type& OperatorDomainMap() const override {}

    /// Method not implemented.
    virtual const super::map_Type& OperatorRangeMap() const override {}

private:
    GetPot                                M_datafile;
    EPETRACOMM                            M_comm;
    BM                                    M_matrix;
    BM                                    M_collapsedMatrix;
    shp<LifeV::BlockEpetra_Map>           M_domainMap;
    shp<LifeV::BlockEpetra_Map>           M_rangeMap;
    shp<BlockMaps>                        M_maps;
};

}

#endif // LINEAROPERATOR_HPP
