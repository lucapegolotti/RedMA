//
// Created by micol on 13/03/2021.
//



#ifndef FSIASSEMBLERFE_HPP
#define FSIASSEMBLERFE_HPP


#include <redma/assemblers/finite_element/NavierStokesAssemblerFE.hpp>
#include <redma/solver/time_marching_algorithms/aTimeMarchingAlgorithm.hpp>
#include <redma/solver/time_marching_algorithms/TimeMarchingAlgorithmFactory.hpp>
namespace RedMA
{
    class FSIAssemblerFE : public NavierStokesAssemblerFE
  {
  public:
      void setup() override;

      FSIAssemblerFE(const DataContainer& data, shp<TreeNode> treeNode,
                              std::string stabilizationName = "");
      void addFSIMassMatrix(shp<aMatrix> mat);
      //void addFSIStiffnessMatrix(shp<aVector> sol, shp<aMatrix> mat);
      shp<aMatrix> getMass(const double& time,
                           const shp<aVector>& sol) override;
      //shp<aMatrix> getMassJacobian(const double& time,const shp<aVector>& sol) override;
      void getBoundaryStiffnessMatrix() const;

      shp<aVector> getRightHandSide(const double& time,
                                    const shp<aVector>& sol) override;

      shp<aMatrix> getJacobianRightHandSide(const double& time,
                                            const shp<aVector>& sol) override;

      void addForcingTerm(const shp<aVector> & rhs  ) const;
      shp<BlockMatrix> returnBoundaryStiffness()const;
      void postProcess(const double& time, const shp<aVector>& sol) override;
      void setExporterDisplacement();
      shp<VECTOREPETRA> getBoundaryIndicator();

  protected:
        void                                              computeLameConstants();

        shp<VECTOREPETRA>                                 M_displacementExporter;

        shp<VECTOREPETRA>                                M_boundaryIndicator;

        shp<BlockMatrix>                                  M_BoundaryStiffness;
        double                                            M_lameI;
        double                                            M_lameII;

        shp<aTimeMarchingAlgorithm>                       M_TMAlgorithm;



  };

  }

  #endif // FSIASSEMBLERFE_HPP
