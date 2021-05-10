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

#ifndef PRECONDITIONEROPERATOREP
#define PRECONDITIONEROPERATOREP

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
class PreconditionerOperator : public LifeV::Operators::LinearOperatorAlgebra
{
    typedef LifeV::Operators::LinearOperatorAlgebra                  super;
    typedef BlockVector                                              BV;
    typedef BlockMatrix                                              BM;

public:
    /// Empty constructor.
    PreconditionerOperator();

    /// Method not implemented.
    virtual int SetUseTranspose(bool UseTranspose) override {}

    /// Method not implemented.
    virtual int Apply(const super::vector_Type& X,
                      super::vector_Type& Y) const override {};

    /*! \brief Apply the approximated inverse to a vector.
     *
     * \param X Vector to which the inverse must be applied.
     * \param Y Result.
     * \return Return code; 0 if successful.
     */
    virtual int ApplyInverse(const super::vector_Type& X,
                             super::vector_Type& Y) const override = 0;

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

    /*! \brief Getter for setup time.
     *
     * \return The setup time.
     */
    inline double getSetupTime() const {return M_setupTime;}

protected:
    EPETRACOMM                      M_comm;
    double                          M_setupTime;
};

}

#endif // PRECONDITIONEROPERATOREP
