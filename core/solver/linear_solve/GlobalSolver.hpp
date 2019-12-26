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

#ifndef GLOBALSOLVER_HPP
#define GLOBALSOLVER_HPP

#include <fstream>

#include <GlobalAssembler.hpp>
#include <TreeStructure.hpp>
#include <GeometryParser.hpp>
#include <TreeStructure.hpp>
#include <GlobalSolverOperator.hpp>

#include <AbstractFunctor.hpp>

#include <lifev/core/array/MapVector.hpp>
#include <lifev/core/array/MapEpetra.hpp>
#include <lifev/core/filter/GetPot.hpp>
#include <TimeMarchingAlgorithmsFactory.hpp>

namespace RedMA
{

class GlobalSolver
{
    typedef std::shared_ptr<Epetra_Comm>               commPtr_Type;
    typedef LifeV::MapEpetra                           map_Type;
	typedef std::shared_ptr<map_Type>                  mapPtr_Type;
    typedef LifeV::MapVector<map_Type>                 MapVector;
    typedef std::shared_ptr<MapVector>                 MapVectorPtr;
    typedef TimeMarchingAlgorithm                      TimeMarchingAlgorithmType;
    typedef std::shared_ptr<TimeMarchingAlgorithmType> TimeMarchingAlgorithmPtr;

    typedef std::function<double(double const&,
                                 double const&,
                                 double const&,
                                 double const&,
                                 unsigned int const& )> FunctionType;

public:
    GlobalSolver(const GetPot& datafile, commPtr_Type comm,
                 bool verbose = false, AbstractFunctor* exactSolution = nullptr);

    void solve();

    void setLawInflow(std::function<double(double)> maxLaw);

    void setLawDtInflow(std::function<double(double)> maxLawDt);

    void setExportNorms(std::string filename);

    void setExportErrors(std::string filename);

    void printMeshSize(std::string filename);

    void setForcingFunction(FunctionType forcingFunction,
                            FunctionType forcingFunctionDt);

private:
    // we pass dt as reference to allow for time adaptvity
    void solveTimestep(const double& time, double& dt);

    void setExactSolution(AbstractFunctor* exactSolution);

    GeometryParser                 M_geometryParser;
    TreeStructure                  M_tree;
    GetPot                         M_datafile;
    commPtr_Type                   M_comm;
    bool                           M_verbose;
    bool                           M_exportNorms;
    bool                           M_exportErrors;
    std::string                    M_filename;
    GlobalAssembler                M_globalAssembler;
    TimeMarchingAlgorithmPtr       M_timeMarchingAlgorithm;
};

}  // namespace RedMA

#endif  // GLOBALSOLVER_HPP
