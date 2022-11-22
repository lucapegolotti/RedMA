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

#ifndef MATRICESGENERATORFIXEDGEOMETRY_HPP
#define MATRICESGENERATORFIXEDGEOMETRY_HPP


#include "aMatricesGenerator.hpp"

namespace RedMA
{

/*! \brief Class for generating the finite element matrices necessary to the
* offline phase of the reduced basis method, for the single block and fixed geometry case
*/
class MatricesGeneratorFixedGeometry : public aMatricesGenerator
{

typedef aAssembler                                      AssemblerType;

public:
    /*! \brief Constructor.
     *
     * \param data A DataContainer.
     * \param comm The MPI Communicator.
     */
    MatricesGeneratorFixedGeometry(const DataContainer& data,
                                   EPETRACOMM comm);

    /*! \brief Generate matrices related to problem (the one specified in the
     * geometry file).
     */
    virtual void generate() override;

protected:
    virtual void createAssemblers() override;

    shp<TreeNode> generateTreeNode();

    shp<AssemblerType>                                  M_assembler;
};

}  // namespace RedMA


#endif //MATRICESGENERATORFIXEDGEOMETRY_HPP
