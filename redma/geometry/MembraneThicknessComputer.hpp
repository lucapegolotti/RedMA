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
#include <redma/problem/DataContainer.hpp>

#include <lifev/core/fem/BCManage.hpp>
#include <lifev/eta/expression/Integrate.hpp>
#include <lifev/core/algorithm/LinearSolver.hpp>
#include <lifev/core/algorithm/PreconditionerML.hpp>
#include <Teuchos_ParameterList.hpp>
#include <Teuchos_XMLParameterListHelpers.hpp>
#include <Teuchos_RCP.hpp>

namespace RedMA {

/*! \brief Class that allows to compute the membrane thickness of the vessel wall in a building block.
 *
 * In the Coupled Momentum method, the structure of the vessel wall is approximated with a thin membrane
 * of a certain thickness. We refer to the following work for further information:
 *
 * Colciago, C.M., Deparis, S. and Quarteroni, A., 2014.
 * Comparisons between reduced order models and full 3D models for fluid–structure interaction
 * problems in haemodynamics.
 * Journal of Computational and Applied Mathematics, 265, pp.120-138.
 *
 * This class allows to compute the membrane thickness over all the vessel wall. This computation can be
 * done in two ways, depending on the value of a suitable flag, called "constant_thickness", set in the datafile.
 * If the flag is true, the thickness equals a fixed percentage (set in the datafile) of the average between
 * the inlet diameter and the average outlet diameters. I the flag is false, the membrane thickness is
 * computed by solving a homogeneous Laplacian problem, where Dirichlet BCs, prescribing a thickness equal
 * to a fixed percentage (set in the datafile) of the vessel diameter, are imposed at the I/O rings.
 *
 */
class MembraneThicknessComputer {

    public:

        /*! \brief Constructor, taking as input a datafile and an Epetra Communicator
         *
         * The constructor initializes the wall flag, the percentage of wall diameter to be used
         * to compute the membrane thickness and a flag telling whether a constant thickness is
         * desired or not. Additionally, a linear system solver is instantiated.
         *
         * @param datafile The datafile
         * @param comm The Epetra Communicator (for multi-processing)
         */
        MembraneThicknessComputer(const DataContainer& datafile,
                                  EPETRACOMM comm);

        /*! \brief Setup method, initializing all quantities needed to compute the membrane thickness.
         *
         * In particular the setup method defines the FE spaces (normal and ETA) for the thickness, a
         * VectorEpetra (set to 0) storing the thickness. Then, if M_constantFlag is true, it just defines
         * the reference radius as the average between the inlet one and the average of the outlet ones.
         * Conversely, if M_constantFlag is false, it assembles the boundary stiffness matrix and the
         * right-hand side term, applying the proper BCs at the I/O rings, and it initializes the linear
         * system solver.
         *
         * @param R_in Vector storing the radia of the inlets
         * @param R_out Vector storing the radia of the outlets
         * @param inletFlag Vector storing the flags of the inlet rings
         * @param outletFlags Vector storing the flags of the outlet rings
         * @param wallFlag Flag of the lateral wall (i.e. structure)
         * @param mesh Shared pointer to the mesh
         */
        void setup(const std::vector<double>& R_in, const std::vector<double>& R_out,
                   const std::vector<unsigned int>& inletFlags, const std::vector<unsigned int>& outletFlags,
                   const unsigned int& wallFlag, shp<MESH> mesh);

        /*! \brief Method to compute the membrane thickness.
         *
         * Depending on the value of M_constantFlag, either a constant thickness is set at the vessel
         * wall or a homogeneous Laplacian problem with suitable BCs at the I/O rings is solved.
         *
         */
        void solve();

        /*! Getter method for the membrane thickness
         *
         * @return Shared pointer to an EpetraVector, storing the membrane thickness of the building block
         */
        shp<VECTOREPETRA> getThickness() const;

    private:

        /*! \brief Method to assemble the stiffness matrix at the vessel wall
         *
         */
        void computeStiffness();

        /*! \brief Method to assemble the right-hand side term for the homogeneous laplacian
         * problem at the vessel wall
         *
         */
        void computeRhs();

        /*! \brief Method to apply the Dirichlet BCs at the I/O rings
         *
         * If M_constantFlag is false, then the membrane thickness is derived by solving a homogeneous
         * laplacian problem. This method imposes the proper Dirichlet BCs at the I/O rings, where
         * the thickness is imposed as equal to a fixed percentage (set at datafile) of the vessel diameter.
         * The default value of such percentage is 10%.
         *
         */
        void applyBCs();

        /*! \brief Method to setup the linear system solver
         *
         */
        void setupLinearSolver();

        /*! \brief Static method, defining a constant (in space and time) function
         *
         * @param t Time
         * @param x x-coordinate
         * @param y y-coordinate
         * @param z z-coordinate
         * @param i Index of the component (0->velocity (vectorial); 1->pressure (scalar))
         * @param K Constant value
         * @return Constant value K, whichever the other inputs
         */
        static double constantFunction(const double &t, const double &x, const double &y,
                                       const double &z, const unsigned int &i,
                                       const double &K);

        shp <VECTOREPETRA>                            M_thickness;
        double                                        M_thicknessProportion;
        int                                           M_constantFlag;

        shp <MATRIXEPETRA>                            M_stiffness;
        shp <VECTOREPETRA>                            M_rhs;

        DataContainer                                 M_datafile;

        EPETRACOMM                                    M_comm;

        LifeV::LinearSolver                           M_linearSolver;
        shp<LifeV::Preconditioner>                    M_preconditioner;

        shp<FESPACE>                                  M_fespace;
        shp<ETFESPACE1>                               M_fespaceETA;

        std::vector<double>                           M_R_in;
        std::vector<double>                           M_R_out;
        double                                        M_R_ref;

        std::vector<unsigned int>                     M_inletFlags;
        std::vector<unsigned int>                     M_outletFlags;
        unsigned int                                  M_wallFlag;

    };

};


#endif //REDMA_MEMBRANETHICKNESSCOMPUTER_H
