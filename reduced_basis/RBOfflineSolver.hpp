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

#ifndef RBOFFLINESOLVER_HPP
#define RBOFFLINESOLVER_HPP

#include <GlobalSolver.hpp>

#include <rb/reduced_basis/fem/ApproximatedAffineFemProblem.hpp>
#include <rb/reduced_basis/rbSolver/ParameterHandler.hpp>

namespace RedMA
{

template<class AssemblerType>
class RBOfflineSolver : public GlobalSolver<AssemblerType>,
                        public rbLifeV::ApproximatedAffineFemProblem
{
    typedef std::shared_ptr<Epetra_Comm>                        commPtr_Type;
    typedef LifeV::RegionMesh<LifeV::LinearTetra>               Mesh;
    typedef std::shared_ptr<Mesh>                               MeshPtr;
    typedef LifeV::FESpace<Mesh, LifeV::MapEpetra>              FESpace;
    typedef std::shared_ptr<FESpace>                            FESpacePtr;
    typedef FE_VECTOR                                           VectorFE;
    typedef std::shared_ptr<VectorFE>                           VectorFEPtr;
    typedef FE_MATRIX                                           MatrixFE;
    typedef std::shared_ptr<MatrixFE>                           MatrixFEPtr;

public:
    RBOfflineSolver(const GetPot& datafile, rbLifeV::ParameterHandler & parameterHandler,
                    commPtr_Type comm, bool verbose = false);

    virtual int finalizeProblem(){};

    virtual void initializeSpecificFemProblem(){};

    virtual void setYh(){};

    virtual const std::shared_ptr<LifeV::MapEpetra> getFieldMap(int iFields) const {};

    virtual FESpacePtr getFieldFeSpace(int _field) const {};

    virtual void setParameter(const param_Type& _mu) {};

    virtual int compute_affine_matrix(const rbLifeV::AffineDecompositionOperator& Op,
                                      const uint& q, FE_MATRIX * & A,
                                      const int _iField = -1,
                                      const int _jField = -1 ) {};

    virtual int compute_affine_vector(const rbLifeV::AffineDecompositionOperator& Op,
                                      const uint& q, VectorFE* & A,
                                      const int _iField = -1 ) {};

    virtual LifeV::MapEpetra const&  globalMap() const {};

    virtual const MeshPtr mesh( ) const {};

    virtual void postprocess(const VectorFE& solution, LifeV::Real time ) {};

    virtual void postprocess(double time = 0.0) {};

    virtual void postprocess(VectorFEPtr * const& solution,
                             LifeV::Real time) {};

    virtual MatrixFEPtr getYh( int iFields ) const {};

    virtual MatrixFEPtr getFieldMatrix(int _iField, int _jField) {};

    virtual void addExporterVariable(std::shared_ptr<LifeV::ExporterHDF5<Mesh> >& _exporter,
                                     bool _addForSupremizer ) {};

    virtual void setExportedSnapshot(VectorFEPtr const * const _exportedVariables ) {};

    virtual void setImportedSnapshot(VectorFEPtr const * const _importedVariables ) {};

    virtual void saveBasisFunctions(const VectorFE& solution,
                                    const LifeV::Real& time) {};

    virtual void initExporter() {};

    virtual void closeExporter() {};

    virtual MatrixFEPtr getMassMatrix(int _field = 0) {};

    virtual void assembleStiffnessMatrix(MatrixFE& _A, const param_Type& _mu, int _numVolumes = 0,
                                         LifeV::UInt * _volumes = nullptr,
                                         bool _integrateSubMesh = false ) {};

    virtual double computeResidual(VectorFEPtr * const _uh, const param_Type& _mu,
                                   VectorFEPtr * const _rh = nullptr,
                                   int _numVolumes = 0, LifeV::UInt * _volumes = nullptr,
                                   bool _integrateSubMesh = false ) {};

    virtual double computeNaturalResidual(VectorFEPtr * const _uh,
                                          const param_Type& _mu) {};

    virtual double computeFieldResidual(int _feField, VectorFEPtr * _uh, const param_Type& _mu,
                                        VectorFEPtr * _rh, int _numVolumes = 0,
                                        LifeV::UInt * _volumes = nullptr,
                                        bool _integrateSubMesh = false ) {};

    virtual double computeError(VectorFEPtr * const _uh, const param_Type& _mu ) {};

    virtual double computeAx(const VectorFEPtr& _uh, const param_Type& _mu,
                             VectorFEPtr& _approxFh, int _iField, int _jField,
                             bool _transpose = false ) {};

    virtual double computeAx(VectorFEPtr * const _uh, const param_Type& _mu,
                             VectorFEPtr * const _approxFh, bool _transpose = false ) {};

    virtual double computeAx(VectorFEPtr * const _uh, VectorFEPtr * const _approxFh,
                             bool _transpose = false ) {};

    virtual int solveFeSystem(VectorFEPtr * const _uh, const param_Type & _mu,
                              double _tol = 1.e-11 ) {};

    virtual int solveFeSystem(VectorFEPtr * const _uh, const param_Type & _mu,
                              VectorFEPtr * const _fh, preconditionerOperatorPtr_Type,
                              double _tol = 1.e-11, int _maxIterations = 100 ) {};

    virtual double assembleRhsPtr(VectorFEPtr * const _fh,  const param_Type & _mu,
                                  int _numVolumes = 0, LifeV::UInt * _volumes = nullptr,
                                  bool _integrateSubMesh = false ) {};

    virtual void setPreconditioner(const param_Type & _mu ) {};

    virtual preconditionerOperatorPtr_Type getPreconditionerOperator() {};

    virtual void applyPreconditioner(VectorFEPtr * const _rh, VectorFEPtr * const _Prh,
                                     const param_Type & _mu) {};

    virtual void applyPreconditioner(VectorFEPtr * const _rh, VectorFEPtr * const _Prh) {};

    virtual linearOperatorPtr_Type getLinearOperator( ) const {};

    virtual void initOfflineExporter(std::string _location) {};

    virtual void closeOfflineExporter() {};

    virtual void applyBCvector(VectorFEPtr& _uh, int _field) {};

    virtual void saveFieldBasisFunctions(VectorFEPtr& _uh, int _field,
                                         const LifeV::Real& time) {};

    virtual void setCollectSnapshot(bool _collectRbSnapshot) {};

    virtual void setUseRbApproximation(bool _useRbApproximation) {};

    virtual LifeV::Real normYh(VectorFEPtr const _u,
                               unsigned int const & _feField) const {};

    virtual double assembleRhsPtr(VectorFEPtr * const _fh,
                                  const param_Type & _mu,
                                  std::vector<VectorFEPtr*> const & _uh,
                                  int _timestep ) {};

    virtual int buildParabolicReducedOperator(const RBmatrix& _An,
                                              const RBmatrix& _Mn,
                                              const param_Type & _mu,
                                              RBmatrix& _timeDependentMatrix,
                                              double _time = 0.0,
                                              RBmatrix * _Rn = nullptr,
                                              RBmatrix * _Nn = nullptr ) {};

    virtual int buildParabolicReducedRhs(RBvector & _timeDependentRhs,
                                         RBvector const & _FN,
                                         RBmatrix const & _AN,
                                         RBmatrix const & _MN,
                                         RBvector const & _mFN,
                                         std::vector<RBvector> const & _UN,
                                         LifeV::Real const & timeStep,
                                         LifeV::Real const & time,
                                         RBvector * _nonlinearLifting = nullptr) {};

    virtual int buildParabolicLsReducedOperator(const RBmatrix& _AAn,
                                                const RBmatrix& _AMn,
                                                const RBmatrix& _MAn,
                                                const RBmatrix& _MMn,
                                                const param_Type & _mu,
                                                RBmatrix& _timeDependentMatrix) {};

    virtual int buildParabolicLsReducedRhs(RBvector & _timeDependentRhs,
                                           RBvector const & _A_F,
                                           RBmatrix const & _A_M,
                                           RBmatrix const & _M_M,
                                           RBvector const & _M_F,
                                           RBvector const & _A_mF,
                                           RBvector const & _M_mF,
                                           RBvector  const & _UN,
                                           LifeV::Real const & time) {};


protected:

};

}  // namespace RedMA

#include "RBOfflineSolver_imp.hpp"

#endif  // RBNAVIERSTOKES_HPP
