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
            loadSTDVector(infile, numMyEntries);
            loadSTDVector(infile, partialSumMyEntries);
            loadSTDVector(infile, myRowMatrixEntriesOfMagicPoints);
            loadSTDVector(infile, myColMatrixEntriesOfMagicPoints);
            loadSTDVector(infile, rowLocalReducedIndices);
            loadSTDVector(infile, globalReducedNodes);
            loadSTDVector(infile, myGlobalReducedNodes);
            loadSTDVector(infile, reducedElements);
            loadSTDVector(infile, numMyReducedEntries);
            loadSTDVector(infile, columnLocalReducedIndices);
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

    void loadMatrix(std::ifstream& infile, std::vector<std::vector<int>>& matrix,
                    int dim1, std::vector<int> dim2)
    {
        matrix.resize(dim1);

        std::string value;
        std::string line;
        getline(infile,line);

        std::stringstream linestream(line);
        getline(linestream,value,','); // we discard the first value (field name)

        for (int i = 0; i < dim1; i++)
        {
            matrix[i].resize(dim2[i]);
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

    template <typename Type>
    void loadSTDVector(std::ifstream& infile, std::vector<Type>& targetVector)
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

        outfile << array2string("numMyEntries", numMyEntries.data(), numMyRows);
        outfile << array2string("partialSumMyEntries", partialSumMyEntries.data(), numMyRows + 1);
        outfile << array2string("myRowMatrixEntriesOfMagicPoints", myRowMatrixEntriesOfMagicPoints.data(), numMyLocalMagicPoints);
        outfile << array2string("myColMatrixEntriesOfMagicPoints", myColMatrixEntriesOfMagicPoints.data(), numMyLocalMagicPoints);
        outfile << array2string("rowLocalReducedIndices", rowLocalReducedIndices.data(), numMyLocalMagicPoints);
        outfile << array2string("globalReducedNodes", globalReducedNodes.data(), 2 * numMyLocalMagicPoints);
        outfile << array2string("myGlobalReducedNodes", myGlobalReducedNodes.data(), numMyGlobalReducedNodes);
        outfile << array2string("reducedElements", reducedElements.data(), numReducedElements);
        outfile << array2string("numMyReducedEntries", numMyReducedEntries.data(), numMyLocalMagicPoints);
        outfile << array2string("columnLocalReducedIndices", columnLocalReducedIndices.data(), numMyLocalMagicPoints);
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
        std::ostringstream streamObj;
        streamObj << value;
        return name + "," + streamObj.str() + "\n";
    }

    template <typename Type>
    std::string array2string(std::string name, Type* array, int numEntries)
    {
        std::ostringstream streamObj;

        std::string retString = name + ",";
        for (int i = 0; i < numEntries; i++)
        {
            streamObj.str() = "";
            streamObj << array[i];
            retString += streamObj.str();
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
        std::ostringstream streamObj;

        std::string retString = name + ",";
        for (int i = 0; i < dim1; i++)
        {
            for (int j = 0; j < dim2; j++)
            {
                streamObj.str() = "";
                streamObj << mat[i][j];
                retString += streamObj.str();
                retString += ",";
            }
        }
        retString = retString.substr(0, retString.size()-1);
        retString += "\n";
        return retString;
    }

    std::string matrix2string(std::string name, std::vector<std::vector<int>> mat,
                              int dim1, std::vector<int> dim2)
    {
        std::ostringstream streamObj;

        std::string retString = name + ",";
        for (int i = 0; i < dim1; i++)
        {
            for (int j = 0; j < dim2[i]; j++)
            {
                streamObj.str() = "";
                streamObj << mat[i][j];
                retString += streamObj.str();
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
    std::vector<int>                numMyEntries;
    std::vector<int>                partialSumMyEntries;
    std::vector<int>                myRowMatrixEntriesOfMagicPoints;
    std::vector<int>                myColMatrixEntriesOfMagicPoints;
    std::vector<int>                rowLocalReducedIndices;
    std::vector<int>                globalReducedNodes;
    std::vector<int>                myGlobalReducedNodes;
    std::vector<unsigned int>       reducedElements;
    std::vector<int>                numMyReducedEntries;
    std::vector<int>                columnLocalReducedIndices;
    std::vector<int>                myLocalMagicPoints;
    std::vector<int>                localIndicesMagicPoints;
    std::vector<int>                magicPointsProcOwner;
    std::vector<int>                globalIndicesMagicPoints;
    std::vector<std::vector<int>>   columnIndices;
    Epetra_SerialDenseMatrix        Qj;
    SHP(MAPEPETRA)                  vectorMap;
};

typedef boost::numeric::ublas::matrix<SHP(MDEIMStructure)>   BlockMDEIMStructure;

}  // namespace RedMA

#endif  // SINGLEMDEIMSTRUCTURE_HPP
