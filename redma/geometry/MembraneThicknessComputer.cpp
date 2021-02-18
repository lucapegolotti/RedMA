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

#include "MembraneThicknessComputer.hpp"

namespace RedMA {

    MembraneThicknessComputer::
    MembraneThicknessComputer(const GetPot& datafile, EPETRACOMM comm):
            M_datafile(datafile), M_comm(comm), M_R_in(1.0)
    {
        M_wallFlag = datafile("structure/flag", 10);
        M_constantFlag = datafile("structure/constant_thickness", 0);
        M_thicknessProportion = datafile("structure/thickness", 0.1);
        M_linearSolver = LifeV::LinearSolver(comm);
    }

    void
    MembraneThicknessComputer::
    setup(const double& R_in, const std::vector<double>& R_out,
          const int& inletFlag, const std::vector<int>& outletFlags,
          shp<MESH> mesh)
    {
        M_fespace.reset(new FESPACE(mesh, "P1", 1, M_comm));
        M_fespaceETA.reset(new ETFESPACE1(M_fespace->mesh(),
                                          &(M_fespace->refFE()),
                                          M_comm));

        M_thickness.reset(new VECTOREPETRA(M_fespace->map(), LifeV::Unique));
        M_thickness->zero();

        if (!M_constantFlag)
        {
            M_R_in = R_in;
            M_R_out = R_out;

            M_inletFlag = inletFlag;
            M_outletFlags = outletFlags;

            this->computeStiffness();
            this->computeRhs();
            this->applyBCs();

            this->setupLinearSolver();
        }
        else
        {
            // If the thickness is chosen as constant, then the reference radius is the average
            // between the inlet radius and the average of the outlet radia.
            double R_out_mean = (std::accumulate(std::begin(R_out), std::end(R_out), 0.0)) / (R_out.size());
            M_R_in = 0.5 * (R_in + R_out_mean);
        }
    }

    void
    MembraneThicknessComputer::
    computeStiffness() {
        using namespace LifeV::ExpressionAssembly;

        M_stiffness.reset(new MATRIXEPETRA(M_fespace->map()));

        LifeV::QuadratureBoundary myBDQR(LifeV::buildTetraBDQR(LifeV::quadRuleTria4pt));

        integrate(boundary(M_fespaceETA->mesh(), M_wallFlag),
                  myBDQR,
                  M_fespaceETA,
                  M_fespaceETA,
                  dot(grad(phi_i), grad(phi_j))
        ) >> M_stiffness;
        M_stiffness->globalAssemble();
    }

    void
    MembraneThicknessComputer::
    computeRhs() {
        M_rhs.reset(new VECTOREPETRA(M_fespace->map(), LifeV::Unique));
        M_rhs->zero();
    }

    void
    MembraneThicknessComputer::
    applyBCs() {
        shp<LifeV::BCHandler> bcs;
        bcs.reset(new LifeV::BCHandler);

        // Imposing the BCs at the rings and at inlet/outlet
        auto BC_in = std::bind(this->constantFunction,
                               std::placeholders::_1,
                               std::placeholders::_2,
                               std::placeholders::_3,
                               std::placeholders::_4,
                               std::placeholders::_5,
                               2.0 * M_R_in * M_thicknessProportion
        );
        LifeV::BCFunctionBase inflowBC(BC_in);

        bcs->addBC("InletRing", M_inletFlag,
                   LifeV::EssentialEdges, LifeV::Full, inflowBC, 3);

        unsigned int count = 0;
        for (double R_out : M_R_out) {
            // TODO: if multiple outlets with different diameters are present, the code must be updated!
            auto BC_out = std::bind(this->constantFunction,
                                    std::placeholders::_1,
                                    std::placeholders::_2,
                                    std::placeholders::_3,
                                    std::placeholders::_4,
                                    std::placeholders::_5,
                                    2.0 * R_out * M_thicknessProportion
            );
            LifeV::BCFunctionBase outflowBC(BC_out);

            bcs->addBC("OutletRing", M_outletFlags[count],
                       LifeV::EssentialEdges, LifeV::Full, outflowBC, 3);

            count++;
        }

        // Applying the defined BCs
        CoutRedirecter ct;
        ct.redirect();

        bcs->bcUpdate(*M_fespace->mesh(), M_fespace->feBd(), M_fespace->dof());
        auto domainMap = M_stiffness->domainMapPtr();
        auto rangeMap = M_stiffness->rangeMapPtr();
        bcManageMatrix(*M_stiffness, *M_fespace->mesh(), M_fespace->dof(),
                       *bcs, M_fespace->feBd(), 1.0, 0.0);
        M_stiffness->globalAssemble(domainMap, rangeMap);

        bcManageRhs(*M_rhs, *M_fespace->mesh(), M_fespace->dof(),
                    *bcs, M_fespace->feBd(), 1.0, 0.0);

        printlog(CYAN, ct.restore(), false);
    }

    void
    MembraneThicknessComputer::
    setupLinearSolver() {

        // Reading solver parameters
        Teuchos::RCP <Teuchos::ParameterList> aztecList = Teuchos::rcp(new Teuchos::ParameterList);
        std::string XMLsolver = this->M_datafile("geometric_structure/xmldeformer",
                                                 "datafiles/SolverParamList.xml");
        aztecList = Teuchos::getParametersFromXmlFile(XMLsolver);
        M_linearSolver.setParameters(*aztecList);

        // Setting up a fake ML preconditioner
        typedef LifeV::PreconditionerML         precML_type;
        precML_type* precRawPtr;
        precRawPtr = new LifeV::PreconditionerML;
        GetPot dummyDatafile;
        precRawPtr->setDataFromGetPot(dummyDatafile, "precMLL");
        M_preconditioner.reset(precRawPtr);
    }

    void
    MembraneThicknessComputer::
    solve() {

        if (!M_constantFlag) {
            // Finalizing linear solver setup
            M_linearSolver.setOperator(M_stiffness);
            M_linearSolver.setPreconditioner(M_preconditioner);
            M_linearSolver.setRightHandSide(M_rhs);

            // Computing the solution
            CoutRedirecter ct;
            ct.redirect();

            M_linearSolver.solve(M_thickness);

            printlog(CYAN, ct.restore(), false);
        }
        else
        {
            shp<VECTOREPETRA> tmp;
            tmp.reset(new VECTOREPETRA(M_fespace->map(), LifeV::Unique));
            tmp->zero();
            std::cout<<"I am here!"<<std::endl;
            *tmp += 1.0;
            *tmp *= (2.0 * M_R_in * M_thicknessProportion);
            *M_thickness += *tmp;
        }
    }

    shp<VECTOREPETRA>
    MembraneThicknessComputer::
    getThickness() const
    {
        return M_thickness;
    }

    double
    MembraneThicknessComputer::
    constantFunction(const double &t, const double &x, const double &y,
                     const double &z, const unsigned int &i, const double &K) {
        return K;
    }

}
