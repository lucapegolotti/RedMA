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

#ifndef MDEIMSTRUCTURE_HPP
#define MDEIMSTRUCTURE_HPP

#include <redma/RedMA.hpp>

namespace RedMA
{

struct MDEIMStructure
{
    MDEIMStructure() :
      allocated(false)
    {

    }

    ~MDEIMStructure()
    {
        if (allocated)
        {
            delete[] numMyEntries;
            delete[] partialSumMyEntries;
            delete[] myRowMatrixEntriesOfMagicPoints;
            delete[] myColMatrixEntriesOfMagicPoints;
            delete[] rowLocalReducedIndices;
            delete[] globalReducedNodes;
            delete[] myGlobalReducedNodes;
            delete[] reducedElements;
            delete[] numMyReducedEntries;
            delete[] columnLocalReducedIndices;

            for (unsigned int i = 0; i < numMyRows; i++)
                delete[] columnIndices[i];

            delete[] columnIndices;
        }
    }

    void dumpMDEIMStructure(std::string dir)
    {
        std::ofstream outfile;
        outfile.open(dir + "/structure.mstr", std::ios_base::out);

        outfile << value2string("N", N);
        outfile << value2string("numGlobalNonzeros", numGlobalNonzeros);
        outfile << value2string("numMyNonzeros", numMyNonzeros);
        outfile << value2string("numMyRows", numMyRows);
        outfile << value2string("numGlobalReducedNodes", numGlobalReducedNodes);
        outfile << value2string("numReducedElements", numReducedElements);
        outfile << value2string("numReducedGlobalNonzeros", numReducedGlobalNonzeros);
        outfile << value2string("numReducedMyNonzeros", numReducedMyNonzeros);
        outfile << value2string("numReducedMyRows", numReducedMyRows);
        outfile << value2string("numMyLocalMagicPoints", numMyLocalMagicPoints);

        outfile << array2string("numMyEntries", numMyEntries, numMyRows);
        outfile << array2string("partialSumMyEntries", partialSumMyEntries, numMyRows + 1);
        outfile << array2string("myRowMatrixEntriesOfMagicPoints", myRowMatrixEntriesOfMagicPoints, numMyLocalMagicPoints);
        outfile << array2string("myColMatrixEntriesOfMagicPoints", myColMatrixEntriesOfMagicPoints, numMyLocalMagicPoints);
        outfile << array2string("rowLocalReducedIndices", rowLocalReducedIndices, numMyLocalMagicPoints);
        outfile << array2string("globalReducedNodes", globalReducedNodes, 2 * numMyLocalMagicPoints);
        outfile << array2string("myGlobalReducedNodes", myGlobalReducedNodes, numMyGlobalReducedNodes);
        outfile << array2string("reducedElements", reducedElements, numReducedElements);
        outfile << array2string("numMyReducedEntries", numMyReducedEntries, numMyLocalMagicPoints);
        outfile << array2string("columnLocalReducedIndices", columnLocalReducedIndices, numMyLocalMagicPoints);
        outfile << array2string("myLocalMagicPoints", myLocalMagicPoints.data(), myLocalMagicPoints.size());
        outfile << array2string("localIndicesMagicPoints", localIndicesMagicPoints.data(), localIndicesMagicPoints.size());
        outfile << array2string("magicPointsProcOwner", magicPointsProcOwner.data(), magicPointsProcOwner.size());
        outfile << array2string("globalIndicesMagicPoints", globalIndicesMagicPoints.data(), globalIndicesMagicPoints.size());

        outfile << matrix2string("Qj", Qj, Qj.RowDim(), Qj.ColDim());
        outfile << matrix2string("columnIndices", columnIndices, numMyRows, numMyEntries);

        outfile.close();
    }

    template <typename Type>
    std::string value2string(std::string name, Type value)
    {
        return name + "," + std::to_string(value) + "\n";
    }

    template <typename Type>
    std::string array2string(std::string name, Type* array, int numEntries)
    {
        std::string retString = name + ",";
        for (int i = 0; i < numEntries; i++)
        {
            retString += std::to_string(array[i]);
            if (i != numEntries - 1)
            {
                retString += ",";
            }
        }
        retString += "\n";
        return retString;
    }

    std::string matrix2string(std::string name, Epetra_SerialDenseMatrix mat,
                              int dim1, int dim2)
    {
        std::string retString = name + ",";
        for (int i = 0; i < dim1; i++)
        {
            for (int j = 0; j < dim2; j++)
            {
                retString += std::to_string(mat[i][j]);
                retString += ",";
            }
        }
        retString = retString.substr(0, retString.size()-1);
        retString += "\n";
        return retString;
    }

    std::string matrix2string(std::string name, int** mat,
                              int dim1, int* dim2)
    {
        std::string retString = name + ",";
        for (int i = 0; i < dim1; i++)
        {
            for (int j = 0; j < dim2[i]; j++)
            {
                retString += std::to_string(mat[i][j]);
                retString += ",";
            }
        }
        retString = retString.substr(0, retString.size()-1);
        retString += "\n";
        return retString;
    }

    int                             N;
    bool                            allocated;
    int                             numGlobalNonzeros;
    int                             numMyNonzeros;
    int                             numMyRows;
    int                             numGlobalReducedNodes;
    int                             numMyGlobalReducedNodes;
    int                             numReducedElements;
    int                             numReducedGlobalNonzeros;
    int                             numReducedMyNonzeros;
    int                             numReducedMyRows;
    int                             numMyLocalMagicPoints;
    int*                            numMyEntries;
    int*                            partialSumMyEntries;
    int*                            myRowMatrixEntriesOfMagicPoints;
    int*                            myColMatrixEntriesOfMagicPoints;
    int*                            rowLocalReducedIndices;
    int*                            globalReducedNodes;
    int*                            myGlobalReducedNodes;
    unsigned int*                   reducedElements;
    int*                            numMyReducedEntries;
    int*                            columnLocalReducedIndices;
    std::vector<int>                myLocalMagicPoints;
    std::vector<int>                localIndicesMagicPoints;
    std::vector<int>                magicPointsProcOwner;
    std::vector<int>                globalIndicesMagicPoints;
    int**                           columnIndices;
    Epetra_SerialDenseMatrix        Qj;
    SHP(MAPEPETRA)                  vectorMap;
};

typedef boost::numeric::ublas::matrix<SHP(MDEIMStructure)>   BlockMDEIMStructure;

}  // namespace RedMA

#endif  // SINGLEMDEIMSTRUCTURE_HPP
