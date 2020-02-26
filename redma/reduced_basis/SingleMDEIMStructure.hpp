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

#ifndef SINGLEMDEIMSTRUCTURE_HPP
#define SINGLEMDEIMSTRUCTURE_HPP

#include <redma/RedMA.hpp>

namespace RedMA
{

struct SingleMDEIMStructure
{
    SingleMDEIMStructure() :
      allocated(false)
    {

    }

    ~SingleMDEIMStructure()
    {
        if (allocated)
        {
            delete[] numMyEntries;
            delete[] partialSumMyEntries;
            delete[] myRowMatrixEntriesOfMagicPoints;
            delete[] myColMatrixEntriesOfMagicPoints;
            delete[] rowLocalReducedIndeces;
            delete[] globalReducedNodes;
            delete[] myGlobalReducedNodes;
            delete[] reducedElements;
            delete[] numMyReducedEntries;
            delete[] columnLocalReducedIndeces;

            for (unsigned int i = 0; i < numMyRows; i++)
                delete[] columnIndeces[i];

            delete[] columnIndeces;
        }
    }

    int                             N;
    bool                            allocated;
    int                             numGlobalNonzeros;
    int                             numMyNonzeros;
    int                             numMyRows;
    int*                            numMyEntries;
    int**                           columnIndeces;
    int*                            partialSumMyEntries;
    int                             numMyLocalMagicPoints;
    std::vector<int>                myLocalMagicPoints;
    std::vector<int>                localIndecesMagicPoints;
    std::vector<int>                magicPointsProcOwner;
    std::vector<int>                globalIndecesMagicPoints;
    int*                            myRowMatrixEntriesOfMagicPoints;
    int*                            myColMatrixEntriesOfMagicPoints;
    int*                            rowLocalReducedIndeces;
    int*                            globalReducedNodes;
    int*                            myGlobalReducedNodes;
    int                             numGlobalReducedNodes;
    int                             numMyGlobalReducedNodes;
    int                             numReducedElements;
    unsigned int*                   reducedElements;
    int                             numReducedGlobalNonzeros;
    int                             numReducedMyNonzeros;
    int                             numReducedMyRows;
    int*                            numMyReducedEntries;
    int*                            columnLocalReducedIndeces;
    Epetra_SerialDenseMatrix        Qj;
    SHP(MAPEPETRA)                  vectorMap;
};

typedef boost::numeric::ublas::matrix<SHP(SingleMDEIMStructure)>   GridMDEIMStructures;

}  // namespace RedMA

#endif  // SINGLEMDEIMSTRUCTURE_HPP
