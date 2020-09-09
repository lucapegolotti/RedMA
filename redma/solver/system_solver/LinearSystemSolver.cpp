#include "LinearSystemSolver.hpp"

namespace RedMA
{

LinearSystemSolver::
LinearSystemSolver(const DataContainer& data) :
  M_data(data),
  M_numSolves(0)
{
}

}
