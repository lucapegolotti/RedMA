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
#include <Exception.hpp>

#include <lifev/eta/fem/ETFESpace.hpp>

#include <functional>

namespace RedMA
{

class NavierStokesAssembler : public AbstractAssembler
{
protected:
    typedef LifeV::ETFESpace<Mesh, MapEpetra, 3, 3>     ETFESpaceVelocity;
    typedef std::shared_ptr<ETFESpaceVelocity>          ETFESpaceVelocityPtr;
    typedef LifeV::ETFESpace<Mesh, MapEpetra, 3, 1>     ETFESpacePressure;
    typedef std::shared_ptr<ETFESpacePressure>          ETFESpacePressurePtr;
    typedef std::function<double(double const&,
                                 double const&,
                                 double const&,
                                 double const&,
                                 unsigned int const& )> FunctionType;


public:
    NavierStokesAssembler(const GetPot& datafile, commPtr_Type comm,
                          const TreeNodePtr& treeNode, bool verbose = false);

    void setup();

    MatrixPtr getMassMatrix(const unsigned int& blockrow,
                            const unsigned int& blockcol);

    MatrixPtr getJacobian(const unsigned int& blockrow,
                          const unsigned int& blockcol);

    inline unsigned int numberOfBlocks() {return 2;}

    void setTimeAndPrevSolution(const double& time,
                                std::vector<VectorPtr> solution);

    std::vector<VectorPtr> computeF();

    std::vector<VectorPtr> computeFder();

protected:
    void assembleConstantMatrices();

    void assembleStiffnessMatrix();

    void assembleDivergenceMatrix();

    void assembleMassMatrix();

    void assembleConvectiveMatrix();

    void assembleJacobianConvectiveMatrix();

    void assembleForcingTerm();

    void assembleForcingTermTimeDerivative();


    FESpacePtr              M_velocityFESpace;
    FESpacePtr              M_pressureFESpace;
    ETFESpaceVelocityPtr    M_velocityFESpaceETA;
    ETFESpacePressurePtr    M_pressureFESpaceETA;

    MatrixPtr               M_A;
    MatrixPtr               M_B;
    MatrixPtr               M_Bt;
    MatrixPtr               M_M;
    MatrixPtr               M_C;
    MatrixPtr               M_J;
    VectorPtr               M_forcingTerm;
    VectorPtr               M_forcingTermTimeDer;
    std::vector<VectorPtr>  M_prevSolution;
    double                  M_time;

    FunctionType            M_forceFunction;
    FunctionType            M_forceTimeDerFunction;
};

}  // namespace RedMA

#endif  // NAVIERSTOKESASSEMBLER_HPP
