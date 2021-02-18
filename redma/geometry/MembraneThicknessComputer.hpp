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

#ifndef REDMA_MEMBRANETHICKNESSCOMPUTER_H
#define REDMA_MEMBRANETHICKNESSCOMPUTER_H

#include <redma/RedMA.hpp>

#include <lifev/core/fem/BCManage.hpp>
#include <lifev/eta/expression/Integrate.hpp>
#include <lifev/core/algorithm/LinearSolver.hpp>
#include <lifev/core/algorithm/PreconditionerML.hpp>
#include <Teuchos_ParameterList.hpp>
#include <Teuchos_XMLParameterListHelpers.hpp>
#include <Teuchos_RCP.hpp>

namespace RedMA {

class MembraneThicknessComputer {

    public:
        MembraneThicknessComputer(const GetPot& datafile,
                                  EPETRACOMM comm);

        void setup(const double& R_in, const std::vector<double>& R_out,
                   const int& inletFlag, const std::vector<int>& outletFlags,
                   shp<MESH> mesh);

        void solve();

        shp<VECTOREPETRA> getThickness() const;

    private:

        void computeStiffness();

        void computeRhs();

        void applyBCs();

        void setupLinearSolver();

        static double constantFunction(const double &t, const double &x, const double &y,
                                       const double &z, const unsigned int &i,
                                       const double &K);

        shp <VECTOREPETRA>                            M_thickness;
        double                                        M_thicknessProportion;
        int                                           M_wallFlag;
        int                                           M_constantFlag;

        shp <MATRIXEPETRA>                            M_stiffness;
        shp <VECTOREPETRA>                            M_rhs;

        GetPot                                        M_datafile;

        EPETRACOMM                                    M_comm;

        LifeV::LinearSolver                           M_linearSolver;
        shp<LifeV::Preconditioner>                    M_preconditioner;

        shp<FESPACE>                                  M_fespace;
        shp<ETFESPACE1>                               M_fespaceETA;

        double                                        M_R_in;
        std::vector<double>                           M_R_out;

        int                                           M_inletFlag;
        std::vector<int>                              M_outletFlags;

    };

};


#endif //REDMA_MEMBRANETHICKNESSCOMPUTER_H
