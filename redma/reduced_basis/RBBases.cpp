#include "RBBases.hpp"

namespace RedMA
{

RBBases::
RBBases(const DataContainer& data, EPETRACOMM comm) :
  M_data(data),
  M_comm(comm),
  M_fespacesAreSet(false)
{
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

    M_onlineTol = 0;
    for (unsigned int i = 0; i < nfields; i++)
    {
        double tolField = M_data("rb/online/basis/podtol_field" + std::to_string(i), 0.0);
        M_onlineTols.push_back(tolField);
        M_onlineTol = M_onlineTol > tolField ? M_onlineTol : tolField;
    }
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
    while (M_onlineTols[index] * M_onlineTols[index] < 1.0 - partialSum / totalenergy)
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
    M_fespacesAreSet = true;
}

void
RBBases::
addVectorsFromFile(std::string filename, std::vector<SHP(VECTOREPETRA)>& vectors,
                   const unsigned int& indexField, int Nmax)
{
    // using namespace boost::filesystem;

    if (std::filesystem::exists(filename))
    {
        std::ifstream infile;
        infile.open(filename);

        unsigned int count = 0;
        std::string line;
        while (std::getline(infile,line) && (Nmax < 0 || count < Nmax))
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
            count++;
        }
        infile.close();
    }
}

void
RBBases::
loadBases()
{
    if (M_fespacesAreSet)
    {
        // load bases
        for (unsigned int i = 0; i < M_numFields; i++)
            addVectorsFromFile(M_path + "/field" + std::to_string(i) + ".basis", M_bases[i], i, M_NsOnline[i]);

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

        for (unsigned int i = 0; i < M_numFields; i++)
        {
            auto curBasis = getEnrichedBasis(i, -1);
            SHP(SparseMatrix) curBasisMatrix(new SparseMatrix(curBasis));
            M_enrichedBasesMatrices.push_back(curBasisMatrix);
            M_enrichedBasesMatricesTransposed.push_back(std::static_pointer_cast<SparseMatrix>(curBasisMatrix->transpose()));
        }
    }
}

void
RBBases::
print()
{
    printlog(WHITE, "\n[RBBases] printing details\n", M_data.getVerbose());
    std::string msg = "";
    if (M_onlineTol > 1e-15)
    {
        for (unsigned int i = 0; i < M_numFields; i++)
            msg += "POD tolerance field " + std::to_string(i) + " = " + std::to_string(M_onlineTols[i]) + "\n";
    }
    else
        msg += "POD tolerance = all vectors\n";
    printlog(WHITE, msg, M_data.getVerbose());

    for (int i = 0; i < M_numFields; i++)
    {
        msg += "Field #" + std::to_string(i) + ":\n";
        msg += "\tLoaded basis size (w/o supremizers) = " + std::to_string(M_bases[i].size()) + "\n";
        if (M_onlineTol > 1e-15)
            msg += "\tSelected basis size = " + std::to_string(M_NsOnline[i]) + "\n";
        for (int j = 0; j < M_numFields; j++)
            msg += "\tPrimal supremizers wrt field " + std::to_string(j) +
                   " = " + std::to_string(M_primalSupremizers(i,j).size()) + "\n";
        msg += "\tDual supremizers = " + std::to_string(M_dualSupremizers[i].size()) + "\n";
        printlog(WHITE, msg, M_data.getVerbose());
        if (M_data("rb/online/basis/useprimalsupremizers", 1))
            msg +=  "\tUsing primal supremizers\n";
        else
            msg +=  "\tNOT using primal supremizers\n";
        printlog(WHITE, msg, M_data.getVerbose());

        if (M_data("rb/online/basis/usedualsupremizers", 0))
        {
            msg +=  "\tUsing dual supremizers\n";
            int nDualSupremizers = M_data("rb/online/basis/numberdualsupremizers", -1);
            if (nDualSupremizers != -1)
                msg +=  "\tUsing " + std::to_string(nDualSupremizers) + " dual supremizers\n";
        }
        else
            msg +=  "\tNOT using dual supremizers\n";
        printlog(WHITE, msg, M_data.getVerbose());

        msg +=  "\tBasis size + supremizers = " + std::to_string(getEnrichedBasis(i,-1).size()) + "\n";
    }
    printlog(WHITE, msg, M_data.getVerbose());
}

void
RBBases::
dump()
{
    // using namespace boost::filesystem;

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

void
RBBases::
scaleBasisWithPiola(unsigned int index, unsigned int ID,
                    std::function<void(SHP(VECTOREPETRA))> transform)
{
    std::vector<SHP(VECTOREPETRA)> newEnrichedBasis;
    auto refEnrichedBasis = getEnrichedBasis(index,-1);

    for (auto vec : refEnrichedBasis)
    {
        SHP(VECTOREPETRA) newVector(new VECTOREPETRA(M_fespaces[index]->map()));
        *newVector = *vec;

        transform(newVector);
        newEnrichedBasis.push_back(newVector);
    }

    M_enrichedBasesMap[index][ID] = newEnrichedBasis;

    SHP(SparseMatrix) basisMatrix(new SparseMatrix(newEnrichedBasis));
    M_enrichedBasesMatricesMap[index][ID] = basisMatrix;
    M_enrichedBasesMatricesTransposedMap[index][ID] = std::static_pointer_cast<SparseMatrix>(basisMatrix->transpose());
}

SHP(BlockMatrix)
RBBases::
leftProject(SHP(BlockMatrix) matrix, unsigned int ID)
{
    SHP(BlockMatrix) projectedMatrix(new BlockMatrix(matrix->nRows(), matrix->nCols()));

    for (unsigned int i = 0; i < matrix->nRows(); i++)
    {
        for (unsigned int j = 0; j < matrix->nCols(); j++)
        {
            projectedMatrix->setBlock(i,j,leftProject(std::static_pointer_cast<SparseMatrix>(matrix->block(i,j)), i, ID));
        }
    }
    return projectedMatrix;
}

SHP(RBMATRIX)
RBBases::
leftProject(SHP(FEMATRIX) matrix, unsigned int basisIndex, unsigned int ID)
{
    SHP(DenseMatrix) retMat(new DenseMatrix);

    if (matrix->data())
    {
        // the basis is transposed

        auto rangeMap = std::static_pointer_cast<MATRIXEPETRA>(getEnrichedBasisMatrices(basisIndex, ID, false)->data())->domainMapPtr();
        auto domainMap = std::static_pointer_cast<MATRIXEPETRA>(matrix->data())->domainMapPtr();

        SHP(MATRIXEPETRA) innerMatrix(new MATRIXEPETRA(*rangeMap));

        // for some reason it is more efficient to transpose and then multiply
        std::static_pointer_cast<MATRIXEPETRA>(getEnrichedBasisMatrices(basisIndex, ID, true)->data())->multiply(false, *std::static_pointer_cast<MATRIXEPETRA>(matrix->data()), false,
                                                                        *innerMatrix);

        innerMatrix->globalAssemble(domainMap, rangeMap);

        SHP(SparseMatrix) retMatEp(new SparseMatrix());
        retMatEp->setMatrix(innerMatrix);

        return retMatEp->toDenseMatrixPtr();
    }
    return retMat;
}

// we assume that the basis for the lagrange multiplier is the identity
SHP(BlockVector)
RBBases::
projectOnLagrangeSpace(SHP(BlockVector) vector)
{
    if (vector->nRows() > 1)
        throw new Exception("projectOnLagrangeSpace: error!");

    SHP(BlockVector) retVec(new BlockVector(1));
    retVec->setBlock(0,std::static_pointer_cast<DistributedVector>(vector->block(0))->toDenseVectorPtr());

    return retVec;
}

SHP(BlockVector)
RBBases::
leftProject(SHP(BlockVector) vector, unsigned int ID)
{
    SHP(BlockVector) projectedVector(new BlockVector(vector->nRows()));

    for (unsigned int i = 0; i < vector->nRows(); i++)
        projectedVector->setBlock(i,leftProject(std::static_pointer_cast<DistributedVector>(vector->block(i)), i, ID));

    return projectedVector;
}

SHP(DenseVector)
RBBases::
leftProject(SHP(DistributedVector) vector, unsigned int basisIndex, unsigned int ID)
{
    if (vector->data())
    {
        SHP(DistributedVector) resVector;
        resVector = std::static_pointer_cast<DistributedVector>(getEnrichedBasisMatrices(basisIndex, ID, true)->multiplyByVector(vector));

        return resVector->toDenseVectorPtr();
    }
    return nullptr;
}

SHP(BlockMatrix)
RBBases::
rightProject(SHP(BlockMatrix) matrix, unsigned int ID)
{
    SHP(BlockMatrix) projectedMatrix(new BlockMatrix(matrix->nRows(), matrix->nCols()));

    for (unsigned int i = 0; i < matrix->nRows(); i++)
    {
        for (unsigned int j = 0; j < matrix->nCols(); j++)
        {
            projectedMatrix->setBlock(i,j,rightProject(std::static_pointer_cast<SparseMatrix>(matrix->block(i,j)), j, ID));
        }
    }
    return projectedMatrix;
}

SHP(DenseMatrix)
RBBases::
rightProject(SHP(FEMATRIX) matrix, unsigned int basisIndex, unsigned int ID)
{
    SHP(DenseMatrix) retMat(new DenseMatrix());

    if (matrix->data())
    {
        // the basis is transposed
        auto rangeMap = std::static_pointer_cast<MATRIXEPETRA>(matrix->data())->rangeMapPtr(); ;
        auto domainMap = std::static_pointer_cast<MATRIXEPETRA>(getEnrichedBasisMatrices(basisIndex, ID, false)->data())->domainMapPtr();

        SHP(MATRIXEPETRA) innerMatrix(new MATRIXEPETRA(*rangeMap));

        std::static_pointer_cast<MATRIXEPETRA>(matrix->data())->multiply(false, *std::static_pointer_cast<MATRIXEPETRA>(getEnrichedBasisMatrices(basisIndex, ID, false)->data()),
                                false, *innerMatrix);

        innerMatrix->globalAssemble(domainMap, rangeMap);

        SHP(SparseMatrix) retMatEp(new SparseMatrix());
        retMatEp->setMatrix(innerMatrix);

        return retMatEp->toDenseMatrixPtr();
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

SHP(aMatrix)
RBBases::
matrixProject(SHP(aMatrix) matrix, unsigned int basisIndexRow,
              unsigned int basisIndexCol,
              unsigned int ID)
{
    SHP(DenseMatrix) retMat(new DenseMatrix());

    if (matrix->data())
    {
        // Vrow' A Vcol

        // compute A Vcol
        auto rangeMap = std::static_pointer_cast<MATRIXEPETRA>(matrix->data())->rangeMapPtr();
        auto domainMap = std::static_pointer_cast<MATRIXEPETRA>(getEnrichedBasisMatrices(basisIndexCol, ID, false)->data())->domainMapPtr();

        SHP(MATRIXEPETRA) auxMatrix(new MATRIXEPETRA(*rangeMap));

        std::static_pointer_cast<MATRIXEPETRA>(matrix->data())->multiply(false, *std::static_pointer_cast<MATRIXEPETRA>(getEnrichedBasisMatrices(basisIndexCol, ID, false)->data()),
                                false, *auxMatrix);
        auxMatrix->globalAssemble(domainMap, rangeMap);

        rangeMap = std::static_pointer_cast<MATRIXEPETRA>(getEnrichedBasisMatrices(basisIndexRow, ID, false)->data())->domainMapPtr();

        SHP(MATRIXEPETRA) innerMatrix(new MATRIXEPETRA(*rangeMap));

        std::static_pointer_cast<MATRIXEPETRA>(getEnrichedBasisMatrices(basisIndexRow, ID, true)->data())->multiply(false, *auxMatrix,
                                                                           false, *innerMatrix);

        innerMatrix->globalAssemble(domainMap, rangeMap);

        SHP(SparseMatrix) retMatEp(new SparseMatrix());
        retMatEp->setMatrix(innerMatrix);

        return retMatEp->toDenseMatrixPtr();
    }
    return retMat;
}

std::vector<SHP(VECTOREPETRA)>
RBBases::
getEnrichedBasis(const unsigned int& index, unsigned int ID)
{
    if (M_enrichedBasesMap[index][ID].size() > 0)
        return M_enrichedBasesMap[index][ID];

    std::vector<SHP(VECTOREPETRA)> retVectors = getBasis(index);
    if (M_data("rb/online/basis/useprimalsupremizers", 1))
    {
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
    }
    if (M_data("rb/online/basis/usedualsupremizers", 0))
    {
        int nDualSupremizers = M_data("rb/online/basis/numberdualsupremizers", -1);
        if (nDualSupremizers == -1)
            for (auto dual : M_dualSupremizers[index])
                retVectors.push_back(dual);
        else
            for (unsigned int i = 0; i < nDualSupremizers; i++)
            {
                if (i < M_dualSupremizers[index].size())
                    retVectors.push_back(M_dualSupremizers[index][i]);
            }
    }
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

SHP(FEMATRIX)
RBBases::
getEnrichedBasisMatrices(const unsigned int& index, const unsigned int& ID,
                         bool transpose)
{
    // if (!transpose)
    // {
    //     if (~M_enrichedBasesMatricesMap[index][ID].isNull())
    //         return M_enrichedBasesMatricesMap[index][ID];
    //
    //     return M_enrichedBasesMatrices[index];
    // }
    // else
    // {
    //     if (~M_enrichedBasesMatricesTransposedMap[index][ID].isNull())
    //         return M_enrichedBasesMatricesTransposedMap[index][ID];
    //
    //     return M_enrichedBasesMatricesTransposed[index];
    // }
    if (!transpose)
    {
        // we have two M_enrichedBasesMatricesMap and M_enrichedBasesMatrices to differentiate
        // between fields that are specific to the domain (e.g. velocity when using piola)
        // and fields that are specific to the building block (e.g. pressure)
        if (M_enrichedBasesMatricesMap[index].find(ID) != M_enrichedBasesMatricesMap[index].end())
            return M_enrichedBasesMatricesMap[index][ID];
        return M_enrichedBasesMatrices[index];
    }
    else
    {
        if (M_enrichedBasesMatricesTransposedMap[index].find(ID) != M_enrichedBasesMatricesTransposedMap[index].end())
            return M_enrichedBasesMatricesTransposedMap[index][ID];
        return M_enrichedBasesMatricesTransposed[index];
    }
}

SHP(VECTOREPETRA)
RBBases::
reconstructFEFunction(SHP(aVector) rbSolution, unsigned int index,
                      unsigned int ID)
{
    SHP(DistributedVector) rbSolutionDistributed;
    if (rbSolution->type() == DENSE)
    {
        rbSolutionDistributed = DistributedVector::convertDenseVector(std::static_pointer_cast<DenseVector>(rbSolution),M_comm);
    }
    else
    {
        rbSolutionDistributed = std::static_pointer_cast<DistributedVector>(rbSolution);
    }


    SHP(DistributedVector) res = std::static_pointer_cast<DistributedVector>(getEnrichedBasisMatrices(index, ID, false)->multiplyByVector(rbSolutionDistributed));

    return std::static_pointer_cast<VECTOREPETRA>(res->data());
}

void
RBBases::
normalizeBasis(const unsigned int& index, SHP(MATRIXEPETRA) normMatrix)
{
    throw new Exception("normalizeBasis not implemented in RBBases");
    // printlog(MAGENTA, "\n[RBBases] normalizing via stabilized Gram Schimdt...\n", M_data.getVerbose());
    //
    // const double thrsh = 1e-5;
    //
    // auto project = [normMatrix] (const SHP(VECTOREPETRA)& vec1,
    //                              const SHP(VECTOREPETRA)& vec2)->double
    // {
    //     auto rmap = *normMatrix->rangeMapPtr();
    //
    //     VECTOREPETRA aux(rmap);
    //     normMatrix->matrixPtr()->Multiply(false, vec1->epetraVector(),
    //                                       aux.epetraVector());
    //
    //     return aux.dot(*vec2);
    // };
    //
    // std::vector<bool> keepVector(M_bases[index].size(), true);
    //
    // auto orthonormalizeWrtBasis = [&](SHP(VECTOREPETRA)& vector,
    //                                   std::vector<SHP(VECTOREPETRA)>& basis,
    //                                   const unsigned int& count)->void
    // {
    //     auto rmap = *normMatrix->rangeMapPtr();
    //
    //     double normVector = project(vector, vector);
    //     *vector /= std::sqrt(normVector);
    //
    //     unsigned int basisIndex = 0;
    //     for (auto& basisV : basis)
    //     {
    //         if (keepVector[basisIndex])
    //         {
    //             // if this is close to 1 then the vectors are almost parallel
    //             // because they are both unitary. Note that we assume basisV to be
    //             // unitary
    //             double coeff = project(vector, basisV);
    //
    //             if (std::abs(1.0 - std::abs(coeff)) > thrsh)
    //             {
    //                 SHP(VECTOREPETRA) aux(new VECTOREPETRA(rmap));
    //                 *aux = *vector - (*basisV) * coeff;
    //
    //                 vector = aux;
    //             }
    //             else
    //             {
    //                 // in this case the supremizer is essentially in the column space.
    //                 // therefore we set it equal to the vector that triggers the check
    //                 std::string msg;
    //                 vector = basisV;
    //                 if (basisIndex > keepVector.size())
    //                 {
    //                     msg = "\nAttention: basisIndex > size of keepVector. ";
    //                     msg += "This indicates that two supremizers are not independent";
    //                     msg += " between each other\n";
    //                     printlog(RED, msg, M_data.getVerbose());
    //                 }
    //                 else
    //                 {
    //                     msg = "\nSwapping supremizer with primal vector\n";
    //                     printlog(WHITE, msg, M_data.getVerbose());
    //                     keepVector[basisIndex] = false;
    //                 }
    //                 return;
    //             }
    //         }
    //         basisIndex++;
    //     }
    //     normVector = project(vector, vector);
    //     *vector /= std::sqrt(normVector);
    // };
    //
    // // re-orthonormalize the primal basis because if it is not orhonormal to
    // // machine precision it's a mess
    //
    // std::vector<SHP(VECTOREPETRA)> incrBasis;
    //
    // printlog(GREEN, "normalizing basis\n", M_data.getVerbose());
    // for (unsigned int i = 0; i < M_bases[index].size(); i++)
    // {
    //     // orthonormalizeWrtBasis(M_bases[index][i], incrBasis, 0);
    //     *M_bases[index][i] /= std::sqrt(project(M_bases[index][i],M_bases[index][i]));
    //     incrBasis.push_back(M_bases[index][i]);
    // }
    //
    // // double check
    // // for (unsigned int i = 10; i < 11; i++)
    // // {
    // //     for (unsigned int j = 0; j < incrBasis.size(); j++)
    // //     {
    // //         double proj = project(incrBasis[i],incrBasis[j]);
    // //         std::cout << "i = " << i << " j = " << j << " proj = " << proj << std::endl;
    // //         if (i == j)
    // //             std::cout << "i = " << i << " j = " << j << " 1 - proj = " << abs(1 - proj) << std::endl;
    // //     }
    // // }
    //
    // // first orthonormalize primal supremizers
    // printlog(GREEN, "normalizing primal supremizers\n", M_data.getVerbose());
    // for (unsigned int j = 0; j < M_numFields; j++)
    // {
    //     if (M_primalSupremizers(index,j).size() > 0)
    //     {
    //         unsigned int count = 0;
    //         for (auto& supr : M_primalSupremizers(index,j))
    //         {
    //             std::string msg = "orthonormalizing primal supremizer ";
    //             msg += std::to_string(count);
    //             msg += "\n";
    //             printlog(YELLOW, msg, M_data.getVerbose());
    //
    //             orthonormalizeWrtBasis(supr, incrBasis, count);
    //             incrBasis.push_back(supr);
    //             count++;
    //         }
    //     }
    // }
    //
    // // then dual supremizers
    // printlog(GREEN, "normalizing dual supremizers\n", M_data.getVerbose());
    // if (M_dualSupremizers[index].size() > 0)
    // {
    //     unsigned int count = 0;
    //     for (auto& supr : M_dualSupremizers[index])
    //     {
    //         std::string msg = "orthonormalizing dual supremizer ";
    //         msg += std::to_string(count);
    //         msg += "\n";
    //         printlog(YELLOW, msg, M_data.getVerbose());
    //
    //         orthonormalizeWrtBasis(supr, incrBasis, count);
    //         incrBasis.push_back(supr);
    //         count++;
    //     }
    // }
    //
    // // finally select the basis functions that we keep
    // std::vector<SHP(VECTOREPETRA)> newBasis;
    // for (unsigned int i = 0; i < keepVector.size(); i++)
    // {
    //     if (keepVector[i])
    //         newBasis.push_back(incrBasis[i]);
    // }
    //
    // M_bases[index] = newBasis;
    //
    // // // double check
    // // for (unsigned int i = 0; i < incrBasis.size(); i++)
    // // {
    // //     for (unsigned int j = 0; j < incrBasis.size(); j++)
    // //     {
    // //         double proj = project(incrBasis[i],incrBasis[j]);
    // //         std::cout << "i = " << i << " j = " << j << " proj = " << proj << std::endl;
    // //     }
    // // }
}

}  // namespace RedMA
