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
    MDEIMStructure(std::string pathfile = "", EPETRACOMM comm = nullptr) :
      allocated(false)
    {
        if (pathfile != "")
        {
            std::ifstream infile;
            infile.open(pathfile, std::ios_base::in);

            // If anyone will ever read this code, I am sorry
            loadValue(infile, N);
            loadValue(infile, numGlobalNonzeros);
            loadValue(infile, numMyNonzeros);
            loadValue(infile, numMyRows);
            loadValue(infile, numGlobalReducedNodes);
            loadValue(infile, numMyGlobalReducedNodes);
            loadValue(infile, numReducedElements);
            loadValue(infile, numReducedGlobalNonzeros);
            loadValue(infile, numReducedMyNonzeros);
            loadValue(infile, numReducedMyRows);
            loadValue(infile, numMyLocalMagicPoints);
            loadValue(infile, Nleft);
            loadValue(infile, Nright);
            loadRawVector(infile, numMyEntries, numMyRows);
            loadRawVector(infile, partialSumMyEntries, numMyRows+1);
            loadRawVector(infile, myRowMatrixEntriesOfMagicPoints, numMyLocalMagicPoints);
            loadRawVector(infile, myColMatrixEntriesOfMagicPoints, numMyLocalMagicPoints);
            loadRawVector(infile, rowLocalReducedIndices, numMyLocalMagicPoints);
            loadRawVector(infile, globalReducedNodes, 2 * numMyLocalMagicPoints);
            loadRawVector(infile, myGlobalReducedNodes, numMyGlobalReducedNodes);
            loadRawVector(infile, reducedElements, numReducedElements);
            loadRawVector(infile, numMyReducedEntries, numMyLocalMagicPoints);
            loadRawVector(infile, columnLocalReducedIndices, numMyLocalMagicPoints);
            loadSTDVector(infile, myLocalMagicPoints);
            loadSTDVector(infile, localIndicesMagicPoints);
            loadSTDVector(infile, magicPointsProcOwner);
            loadSTDVector(infile, globalIndicesMagicPoints);
            loadMatrix(infile, Qj, N, N);
            loadMatrix(infile, columnIndices, numMyRows, numMyEntries);

            int * myGlobalElements = 0;

            vectorMap.reset(new MAPEPETRA(numGlobalNonzeros,
                                          numMyNonzeros,
                                          myGlobalElements, comm));

            infile.close();
            allocated = true;
        }
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

    void loadMatrix(std::ifstream& infile, Epetra_SerialDenseMatrix& matrix,
                    int dim1, int dim2)
    {
        matrix.Reshape(dim1, dim2);

        std::string value;
        std::string line;
        getline(infile,line);

        std::stringstream linestream(line);
        getline(linestream,value,','); // we discard the first value (field name)

        for (int i = 0; i < dim1; i++)
        {
            for (int j = 0; j < dim2; j++)
            {
                getline(linestream, value, ',');
                matrix[i][j] = std::atof(value.c_str());
            }
        }
    }

    void loadMatrix(std::ifstream& infile, int**& matrix, int dim1, int* dim2)
    {
        matrix = new int*[dim1];

        std::string value;
        std::string line;
        getline(infile,line);

        std::stringstream linestream(line);
        getline(linestream,value,','); // we discard the first value (field name)

        for (int i = 0; i < dim1; i++)
        {
            matrix[i] = new int[dim2[i]];
            for (int j = 0; j < dim2[i]; j++)
            {
                getline(linestream, value, ',');
                matrix[i][j] = std::atoi(value.c_str());
            }
        }
    }

    template <typename Type>
    void loadRawVector(std::ifstream& infile, Type*& targetVector, unsigned int numElements)
    {
        targetVector = new Type[numElements];

        std::string line;
        std::getline(infile,line);

        unsigned int count = 0;
        std::stringstream linestream(line);
        std::string value;
        while (getline(linestream,value,','))
        {
            if (count > 0)
                targetVector[count-1] = std::atoi(value.c_str());
            count++;
        }
    }

    void loadSTDVector(std::ifstream& infile, std::vector<int>& targetVector)
    {
        std::string line;
        std::getline(infile,line);

        unsigned int count = 0;
        std::stringstream linestream(line);
        std::string value;
        while (getline(linestream,value,','))
        {
            if (count > 0)
                targetVector.push_back(std::atoi(value.c_str()));
            count++;
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
        outfile << value2string("numMyGlobalReducedNodes", numMyGlobalReducedNodes);
        outfile << value2string("numReducedElements", numReducedElements);
        outfile << value2string("numReducedGlobalNonzeros", numReducedGlobalNonzeros);
        outfile << value2string("numReducedMyNonzeros", numReducedMyNonzeros);
        outfile << value2string("numReducedMyRows", numReducedMyRows);
        outfile << value2string("numMyLocalMagicPoints", numMyLocalMagicPoints);
        outfile << value2string("Nleft", Nleft);
        outfile << value2string("Nright", Nright);

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
    void loadValue(std::ifstream& infile, Type& value)
    {
        std::string line;
        std::getline(infile,line);

        unsigned int count = 0;
        std::stringstream linestream(line);
        std::string stringValue;
        while (getline(linestream,stringValue,','))
        {
            if (count == 1)
                value = std::atoi(stringValue.c_str());
            count++;
        }
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
    int                             Nleft;
    int                             Nright;
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