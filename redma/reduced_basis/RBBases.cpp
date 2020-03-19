#include "RBBases.hpp"

namespace RedMA
{

RBBases::
RBBases(const DataContainer& data, EPETRACOMM comm) :
  M_data(data),
  M_comm(comm)
{
    M_onlineTol = M_data("rb/online/basis/podtol", 0.0);
}

void
RBBases::
setNumberOfFields(const unsigned int& nfields)
{
    M_numFields = nfields;
    M_bases.resize(nfields);
    M_svs.resize(nfields);
    M_NsOnline.resize(nfields);
    M_fespaces.resize(nfields);
    M_primalSupremizers.resize(nfields,nfields);
    M_dualSupremizers.resize(nfields);
}

void
RBBases::
setPath(std::string path)
{
    M_path = path;
}

void
RBBases::
loadSingularValues()
{
    for (unsigned int i = 0; i < M_numFields; i++)
    {
        std::ifstream infile;

        infile.open(M_path + "/svd" + std::to_string(i) + ".txt");

        std::string line;
        while (std::getline(infile, line))
        {
            double newValue = std::atof(line.c_str());
            M_svs[i].push_back(newValue);
        }

        infile.close();

        if (M_onlineTol > 1e-15)
            computeOnlineNumberBasisFunctions(i);
    }
}

void
RBBases::
addPrimalSupremizer(SHP(VECTOREPETRA) supremizer,
                    const unsigned int& fieldToAugment,
                    const unsigned int& fieldConstraint)
{
    M_primalSupremizers(fieldToAugment,fieldConstraint).push_back(supremizer);
}

void
RBBases::
computeOnlineNumberBasisFunctions(unsigned int index)
{
    double totalenergy = 0;

    for (auto sv : M_svs[index])
        totalenergy += sv * sv;

    double partialSum = 0;
    unsigned int N = 0;
    while (M_onlineTol * M_onlineTol < 1.0 - partialSum / totalenergy)
    {
        partialSum += M_svs[index][N] * M_svs[index][N];
        N++;
    }

    M_NsOnline[index] = N;
}

void
RBBases::
addDualSupremizer(SHP(VECTOREPETRA) supremizer, const unsigned int& fieldToAugment)
{
    M_dualSupremizers[fieldToAugment].push_back(supremizer);
}

void
RBBases::
setFESpace(SHP(FESPACE) fespace, const unsigned int& indexbasis)
{
    M_fespaces[indexbasis] = fespace;
}

void
RBBases::
addVectorsFromFile(std::string filename, std::vector<SHP(VECTOREPETRA)>& vectors,
                   const unsigned int& indexField)
{
    using namespace boost::filesystem;

    if (exists(filename))
    {
        std::ifstream infile;
        infile.open(filename);

        std::string line;
        while (std::getline(infile,line))
        {
            SHP(VECTOREPETRA) newVector(new VECTOREPETRA(M_fespaces[indexField]->map()));

            std::stringstream linestream(line);
            std::string value;
            unsigned int index = 0;
            while(getline(linestream,value,','))
            {
                newVector->operator[](index) = std::atof(value.c_str());
                index++;
            }
            if (index != newVector->epetraVector().GlobalLength())
                throw new Exception("Stored basis length does not match fespace dimension!");

            vectors.push_back(newVector);
        }
        infile.close();
    }
}

void
RBBases::
loadBases()
{

    // load bases
    for (unsigned int i = 0; i < M_numFields; i++)
        addVectorsFromFile(M_path + "/field" + std::to_string(i) + ".basis", M_bases[i], i);

    // load primal supremizers
    for (unsigned int i = 0; i < M_numFields; i++)
    {
        for (unsigned int j = 0; j < M_numFields; j++)
        {
            std::string name = "/primal_supremizers_" + std::to_string(i) + "_" + std::to_string(j);
            addVectorsFromFile(M_path + name + ".basis", M_primalSupremizers(i,j), i);
        }
    }

    // load bases
    for (unsigned int i = 0; i < M_numFields; i++)
    {
        addVectorsFromFile(M_path + "/dual_supremizers" + std::to_string(i) + ".basis",
                           M_dualSupremizers[i], i);

    }

    print();

    for (unsigned i = 0; i < M_numFields; i++)
    if (M_bases[i].size() < M_NsOnline[i])
        throw new Exception("Selected tolerance requires more vectors than the ones stored");

}

void
RBBases::
print()
{
    printlog(WHITE, "\n[RBBases] printing details\n", M_data.getVerbose());
    std::string msg = "";
    if (M_onlineTol > 1e-15)
        msg += "POD tolerance = " + std::to_string(M_onlineTol) + "\n";
    else
        msg += "POD tolerance = all vectors\n";

    for (int i = 0; i < M_numFields; i++)
    {
        msg += "Field #" + std::to_string(i) + ":\n";
        msg += "\tTotal basis size = " + std::to_string(M_bases[i].size()) + "\n";
        if (M_onlineTol > 1e-15)
            msg += "\tSelected basis size = " + std::to_string(M_NsOnline[i]) + "\n";
        for (int j = 0; j < M_numFields; j++)
            msg += "\tPrimal supremizers wrt field " + std::to_string(j) +
                   " = " + std::to_string(M_primalSupremizers(i,j).size()) + "\n";
        msg += "\tDual supremizers = " + std::to_string(M_dualSupremizers[i].size()) + "\n";
    }
    printlog(WHITE, msg, M_data.getVerbose());
}

void
RBBases::
dump()
{
    using namespace boost::filesystem;

    bool binary = M_data("rb/basis/dumpbinary", true);

    std::ios_base::openmode omode = std::ios_base::app;
    if (binary)
        omode = omode | std::ios::binary;

    unsigned int fieldIndex = 0;
    for (auto singlebasis : M_bases)
    {
        std::ofstream outfile;
        outfile.open(M_path + "/field" + std::to_string(fieldIndex) + ".basis", omode);

        for (auto vec : singlebasis)
        {
            FEVECTOR curvec;
            curvec.data() = vec;
            std::string str2write = curvec.getString(',') + "\n";

            outfile.write(str2write.c_str(), str2write.size());
        }

        outfile.close();
        fieldIndex++;
    }

    // dump primal supremizers
    for (unsigned int i = 0; i < M_numFields; i++)
    {
        for (unsigned int j = 0; j < M_numFields; j++)
        {
            if (M_primalSupremizers(i,j).size() > 0)
            {
                std::ofstream outfile;
                std::string name = "/primal_supremizers_" + std::to_string(i) + "_" + std::to_string(j);
                outfile.open(M_path + name + ".basis", omode);

                for (auto vec : M_primalSupremizers(i,j))
                {
                    FEVECTOR curvec;
                    curvec.data() = vec;
                    std::string str2write = curvec.getString(',') + "\n";

                    outfile.write(str2write.c_str(), str2write.size());
                }

                outfile.close();
            }
        }
    }

    fieldIndex = 0;
    // dump dual supremizers
    for (auto singlebasis : M_dualSupremizers)
    {
        if (singlebasis.size() > 0)
        {
            std::ofstream outfile;
            outfile.open(M_path + "/dual_supremizers" + std::to_string(fieldIndex) + ".basis", omode);

            for (auto vec : singlebasis)
            {
                FEVECTOR curvec;
                curvec.data() = vec;
                std::string str2write = curvec.getString(',') + "\n";

                outfile.write(str2write.c_str(), str2write.size());
            }
            outfile.close();
        }
        fieldIndex++;
    }
}

BlockMatrix<DenseMatrix>
RBBases::
leftProject(BlockMatrix<MatrixEp> matrix)
{
    BlockMatrix<DenseMatrix> projectedMatrix(matrix.nRows(), matrix.nCols());

    for (unsigned int i = 0; i < matrix.nRows(); i++)
    {
        for (unsigned int j = 0; j < matrix.nCols(); j++)
        {
            projectedMatrix.block(i,j) = leftProject(matrix.block(i,j), i);
        }
    }
    return projectedMatrix;
}

DenseMatrix
RBBases::
leftProject(MatrixEp matrix, unsigned int basisIndex)
{
    DenseMatrix retMat;

    if (matrix.data())
    {
        std::vector<SHP(VECTOREPETRA)> basis = getEnrichedBasis(basisIndex);

        unsigned int nrows = basis.size();
        unsigned int ncols = matrix.data()->domainMapPtr()->mapSize();
        SHP(DENSEMATRIX) innerMatrix(new DENSEMATRIX(nrows, ncols));

        VECTOREPETRA result(*matrix.data()->domainMapPtr());
        for (unsigned int i = 0; i < nrows; i++)
        {
            result.zero();
            matrix.data()->matrixPtr()->Multiply(true, basis[i]->epetraVector(),
                                                 result.epetraVector());

            for (unsigned int j = 0; j < ncols; j++)
                (*innerMatrix)(i,j) = result[j];
        }
        retMat.data() = innerMatrix;
    }
    return retMat;
}

// we assume that the basis for the lagrange multiplier is the identity
BlockVector<DenseVector>
RBBases::
projectOnLagrangeSpace(BlockVector<VectorEp> vector)
{
    if (vector.nRows() > 1)
        throw new Exception("projectOnLagrangeSpace: error!");

    BlockVector<DenseVector> retVec(1);
    retVec.block(0).data() = vector.block(0).toDenseVector().data();

    return retVec;
}

BlockVector<DenseVector>
RBBases::
leftProject(BlockVector<VectorEp> vector)
{
    BlockVector<DenseVector> projectedVector(vector.nRows());

    for (unsigned int i = 0; i < vector.nRows(); i++)
        projectedVector.block(i) = leftProject(vector.block(i), i);

    return projectedVector;
}

DenseVector
RBBases::
leftProject(VectorEp vector, unsigned int basisIndex)
{
    DenseVector retVec;

    if (vector.data())
    {
        std::vector<SHP(VECTOREPETRA)> basis = getEnrichedBasis(basisIndex);

        unsigned int nrows = basis.size();
        SHP(DENSEVECTOR) innerVector(new DENSEVECTOR(nrows));

        double result = 0;
        for (unsigned int i = 0; i < nrows; i++)
        {
            (*innerVector)(i) = basis[i]->dot(*vector.data());
        }
        retVec.data() = innerVector;
    }
    return retVec;
}

BlockMatrix<DenseMatrix>
RBBases::
rightProject(BlockMatrix<MatrixEp> matrix)
{
    BlockMatrix<DenseMatrix> projectedMatrix(matrix.nRows(), matrix.nCols());

    for (unsigned int i = 0; i < matrix.nRows(); i++)
    {
        for (unsigned int j = 0; j < matrix.nCols(); j++)
        {
            projectedMatrix.block(i,j) = rightProject(matrix.block(i,j), j);
        }
    }
    return projectedMatrix;
}

DenseMatrix
RBBases::
rightProject(MatrixEp matrix, unsigned int basisIndex)
{
    DenseMatrix retMat;
    if (matrix.data())
    {
        std::vector<SHP(VECTOREPETRA)> basis = getEnrichedBasis(basisIndex);

        unsigned int nrows = matrix.data()->rangeMapPtr()->mapSize();
        unsigned int ncols = basis.size();
        SHP(DENSEMATRIX) innerMatrix(new DENSEMATRIX(nrows, ncols));

        VECTOREPETRA result(*matrix.data()->rangeMapPtr());

        for (unsigned int j = 0; j < ncols; j++)
        {
            result.zero();
            matrix.data()->matrixPtr()->Multiply(false, basis[j]->epetraVector(),
                                                 result.epetraVector());

            for (unsigned int i = 0; i < nrows; i++)
                (*innerMatrix)(i,j) = result[i];
        }
        retMat.data() = innerMatrix;
    }
    return retMat;
}


std::vector<SHP(VECTOREPETRA)>
RBBases::
getBasis(const unsigned int& index)
{
    if (M_onlineTol < 1e-15)
        return M_bases[index];

    std::vector<SHP(VECTOREPETRA)> basis;

    for (unsigned int i = 0; i < M_NsOnline[index]; i++)
        basis.push_back(M_bases[index][i]);

    return basis;
}

DenseMatrix
RBBases::
matrixProject(MatrixEp matrix, unsigned int basisIndexRow,
                               unsigned int basisIndexCol)
{
    DenseMatrix retMat;
    if (matrix.data())
    {
        std::vector<SHP(VECTOREPETRA)> basisLeft = getEnrichedBasis(basisIndexRow);
        std::vector<SHP(VECTOREPETRA)> basisRight = getEnrichedBasis(basisIndexCol);

        unsigned int Nleft = basisLeft.size();
        unsigned int Nright = basisRight.size();

        SHP(DENSEMATRIX) newMatrix(new DENSEMATRIX(Nleft, Nright));

        SHP(MATRIXEPETRA) matrixEpetra = matrix.data();

        for (unsigned int i = 0; i < Nleft; i++)
        {
            VECTOREPETRA aux(*basisRight[0]->mapPtr());
            matrixEpetra->matrixPtr()->Multiply(true,
                            basisLeft[i]->epetraVector(), aux.epetraVector());

            for (unsigned int j = 0; j < Nright; j++)
            {
                (*newMatrix)(i,j) += aux.dot(*basisRight[j]);
            }
        }

        retMat.data() = newMatrix;
    }
    return retMat;
}

std::vector<SHP(VECTOREPETRA)>
RBBases::
getEnrichedBasis(const unsigned int& index)
{
    std::vector<SHP(VECTOREPETRA)> retVectors = getBasis(index);
    for (unsigned int j = 0; j < M_numFields; j++)
    {
        if (M_onlineTol > 1e-15 && M_primalSupremizers(index,j).size() > 0)
        {
            for (unsigned int i = 0; i < M_NsOnline[j]; i++)
            {
                retVectors.push_back(M_primalSupremizers(index,j)[i]);
            }
        }
        else
        {
            for (auto primal : M_primalSupremizers(index,j))
                retVectors.push_back(primal);
        }
    }

    for (auto dual : M_dualSupremizers[index])
        retVectors.push_back(dual);

    return retVectors;
}

std::vector<unsigned int>
RBBases::
getSelectors(unsigned int index)
{
    std::vector<unsigned int> selectors;

    if (M_onlineTol > 1e-15)
    {
        for (unsigned int i = 0; i < M_NsOnline[index]; i++)
            selectors.push_back(i);

        unsigned int offset = M_bases[index].size();
        // primal supremizers also depend on the ns online
        for (unsigned int j = 0; j < M_numFields; j++)
        {
            if (M_primalSupremizers(index,j).size() > 0)
            {
                for (unsigned int jj = 0; jj < M_NsOnline[j]; jj++)
                    selectors.push_back(jj + offset);
            }
            offset += M_primalSupremizers(index,j).size();
        }

        for (unsigned int i = 0; i < M_dualSupremizers[index].size(); i++)
            selectors.push_back(i + offset);
    }
    return selectors;
}

void
RBBases::
reconstructFEFunction(SHP(VECTOREPETRA) feFunction, DenseVector rbSolution, unsigned int index)
{
    feFunction->zero();

    auto basis = getEnrichedBasis(index);
    for (unsigned int i = 0; i < rbSolution.getNumRows(); i++)
        *feFunction += (*basis[i]) * (*rbSolution.data())(i);
}

void
RBBases::
normalizeBasis(const unsigned int& index, SHP(MATRIXEPETRA) normMatrix)
{
    printlog(MAGENTA, "\n[RBBases] normalizing via stabilized Gram Schimdt...\n", M_data.getVerbose());

    const double thrsh = 1e-5;

    auto project = [normMatrix] (const SHP(VECTOREPETRA)& vec1,
                                 const SHP(VECTOREPETRA)& vec2)->double
    {
        auto rmap = *normMatrix->rangeMapPtr();

        VECTOREPETRA aux(rmap);
        normMatrix->matrixPtr()->Multiply(false, vec1->epetraVector(),
                                          aux.epetraVector());

        return aux.dot(*vec2);
    };

    std::vector<bool> keepVector(M_bases[index].size(), true);

    auto orthonormalizeWrtBasis = [&](SHP(VECTOREPETRA)& vector,
                                      std::vector<SHP(VECTOREPETRA)>& basis,
                                      const unsigned int& count)->void
    {
        auto rmap = *normMatrix->rangeMapPtr();

        double normVector = project(vector, vector);
        *vector /= std::sqrt(normVector);

        unsigned int basisIndex = 0;
        for (auto& basisV : basis)
        {
            if (keepVector[basisIndex])
            {
                // if this is close to 1 then the vectors are almost parallel
                // because they are both unitary
                double coeff = project(vector, basisV) / project(basisV, basisV);

                // std::ostringstream streamOb;
                // streamOb << std::abs(1.0 - std::abs(coeff));
                //
                // std::string msg = std::to_string(count) + ": abs(1 - abs(coeff projection)) = ";
                // msg += streamOb.str();
                // printlog(YELLOW, msg, M_data.getVerbose());

                if (std::abs(1.0 - std::abs(coeff)) > thrsh)
                {
                    SHP(VECTOREPETRA) aux(new VECTOREPETRA(rmap));
                    *aux = *vector - (*basisV) * coeff;
                    double normOrth = sqrt(project(aux, aux));

                    // std::string msg = " norm orthogonal projection = ";
                    // msg += std::to_string(normOrth);
                    // msg += "\n";
                    // printlog(YELLOW, msg, M_data.getVerbose());

                    *aux /= normOrth;
                    vector = aux;
                }
                else
                {
                    // in this case the supremizer is essentially in the column space.
                    // therefore we set it equal to the vector that triggers the check
                    std::string msg;
                    vector = basisV;
                    if (basisIndex > keepVector.size())
                    {
                        msg = "\nAttention: basisIndex > size of keepVector. ";
                        msg += "This indicates that two supremizers are not independent";
                        msg += " between each other\n";
                        printlog(RED, msg, M_data.getVerbose());
                    }
                    else
                    {
                        msg = "\nSwapping supremizer with primal vector\n";
                        printlog(WHITE, msg, M_data.getVerbose());
                        keepVector[basisIndex] = false;
                    }
                    return;
                }
            }
            basisIndex++;
        }
    };

    // re-orthonormalize the primal basis because if it is not orhonormal to
    // machine precision it's a mess

    std::vector<SHP(VECTOREPETRA)> incrBasis;

    printlog(GREEN, "normalizing basis\n", M_data.getVerbose());
    for (unsigned int i = 0; i < M_bases[index].size(); i++)
    {
        orthonormalizeWrtBasis(M_bases[index][i], incrBasis, 0);
        incrBasis.push_back(M_bases[index][i]);
    }

    // double check
    for (unsigned int i = 10; i < 11; i++)
    {
        for (unsigned int j = 0; j < incrBasis.size(); j++)
        {
            double proj = project(incrBasis[i],incrBasis[j]);
            std::cout << "i = " << i << " j = " << j << " proj = " << proj << std::endl;
            if (i == j)
                std::cout << "i = " << i << " j = " << j << " 1 - proj = " << abs(1 - proj) << std::endl;
        }
    }

    // first orthonormalize primal supremizers
    printlog(GREEN, "normalizing primal supremizers\n", M_data.getVerbose());
    for (unsigned int j = 0; j < M_numFields; j++)
    {
        if (M_primalSupremizers(index,j).size() > 0)
        {
            unsigned int count = 0;
            for (auto& supr : M_primalSupremizers(index,j))
            {
                std::string msg = "orthonormalizing primal supremizer ";
                msg += std::to_string(count);
                msg += "\n";
                printlog(YELLOW, msg, M_data.getVerbose());

                orthonormalizeWrtBasis(supr, incrBasis, count);
                incrBasis.push_back(supr);
                count++;
            }
        }
    }

    // then dual supremizers
    printlog(GREEN, "normalizing dual supremizers\n", M_data.getVerbose());
    if (M_dualSupremizers[index].size() > 0)
    {
        unsigned int count = 0;
        for (auto& supr : M_dualSupremizers[index])
        {
            std::string msg = "orthonormalizing primal supremizer ";
            msg += std::to_string(count);
            msg += "\n";
            printlog(YELLOW, msg, M_data.getVerbose());

            orthonormalizeWrtBasis(supr, incrBasis, count);
            incrBasis.push_back(supr);
            count++;
        }
    }

    // finally select the basis functions that we keep
    std::vector<SHP(VECTOREPETRA)> newBasis;
    for (unsigned int i = 0; i < keepVector.size(); i++)
    {
        if (keepVector[i])
            newBasis.push_back(incrBasis[i]);
    }

    M_bases[index] = newBasis;

    // // double check
    // for (unsigned int i = 0; i < incrBasis.size(); i++)
    // {
    //     for (unsigned int j = 0; j < incrBasis.size(); j++)
    //     {
    //         double proj = project(incrBasis[i],incrBasis[j]);
    //         std::cout << "i = " << i << " j = " << j << " proj = " << proj << std::endl;
    //     }
    // }
}

}  // namespace RedMA
