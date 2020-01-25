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

#ifndef LINEAROPERATOREP_HPP
#define LINEAROPERATOREP_HPP

#include <redma/RedMA.hpp>
#include <redma/solver/array/BlockVector.hpp>
#include <redma/solver/array/BlockMatrix.hpp>

#include <lifev/core/filter/GetPot.hpp>
#include <lifev/core/linear_algebra/LinearOperatorAlgebra.hpp>
#include <lifev/core/linear_algebra/BlockEpetra_MultiVector.hpp>
#include <lifev/core/linear_algebra/BlockEpetra_Map.hpp>

namespace RedMA
{

class LinearOperatorEp : public LifeV::Operators::LinearOperatorAlgebra
{
    typedef LifeV::Operators::LinearOperatorAlgebra     super;
    typedef BlockVector<VectorEp>                       BV;
    typedef BlockMatrix<MatrixEp>                       BM;

public:
    LinearOperatorEp();

    void setup(BM matrix, EPETRACOMM comm);

    // I provide null implementation of virtual methods
    // only to be able to instantiate class
    virtual int SetUseTranspose(bool UseTranspose) override {}

    virtual int Apply(const super::vector_Type& X, super::vector_Type& Y) const override;

    virtual int ApplyInverse(const super::vector_Type& X,
                             super::vector_Type& Y) const override {return -1;}

    virtual double NormInf() const override {return -1;}

    virtual const char * Label() const override {}

    virtual bool UseTranspose() const override {return false;}

    virtual bool HasNormInf() const override {return false;}

    virtual const super::comm_Type& Comm() const override {return *M_comm;};

    virtual const super::map_Type& OperatorDomainMap() const override {}

    virtual const super::map_Type& OperatorRangeMap() const override {}

private:
    GetPot                          M_datafile;
    EPETRACOMM                      M_comm;
    BM                              M_matrix;
    SHP(LifeV::BlockEpetra_Map)     M_domainMap;
    SHP(LifeV::BlockEpetra_Map)     M_rangeMap;
};

}

#endif // LINEAROPERATOREP_HPP
