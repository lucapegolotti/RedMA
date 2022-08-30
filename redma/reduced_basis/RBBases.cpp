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
    if (!fs::exists(path))
        throw new Exception("Path " + path + " does not exist! Maybe you forgot to update the .xml geometry file.");

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
        else
            M_NsOnline[i] = std::numeric_limits<unsigned int>::max();
    }
}

void
RBBases::
addPrimalSupremizer(shp<VECTOREPETRA> supremizer,
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
addDualSupremizer(shp<VECTOREPETRA> supremizer, const unsigned int& fieldToAugment)
{
    M_dualSupremizers[fieldToAugment].push_back(supremizer);
}

void
RBBases::
setFESpace(shp<FESPACE> fespace, const unsigned int& indexbasis)
{
    M_fespaces[indexbasis] = fespace;
    M_fespacesAreSet = true;
}

void
RBBases::
addVectorsFromFile(std::string filename, std::vector<shp<VECTOREPETRA>>& vectors,
                   const unsigned int& indexField, int Nmax)
{
    if (fs::exists(filename))
    {
        std::ifstream infile;
        infile.open(filename);

        unsigned int count = 0;
        std::string line;
        while (std::getline(infile,line) && (Nmax < 0 || count < Nmax))
        {
            shp<VECTOREPETRA> newVector(new VECTOREPETRA(M_fespaces[indexField]->map()));

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
            for (unsigned int j = 0; j < M_numFields; j++)
            {
                std::string name = "/primal_supremizers_" + std::to_string(i) + "_" + std::to_string(j);
                addVectorsFromFile(M_path + name + ".basis", M_primalSupremizers(i,j), i);
            }

        // load dual supremizers
        for (unsigned int i = 0; i < M_numFields; i++)
            addVectorsFromFile(M_path + "/dual_supremizers" + std::to_string(i) + ".basis",
                               M_dualSupremizers[i], i);

        /*for (unsigned i = 0; i < M_numFields; i++)
            if (M_bases[i].size() < M_NsOnline[i])
                throw new Exception("Selected tolerance requires more vectors than those stored");*/

        for (unsigned int i = 0; i < M_numFields; i++)
        {
            auto curBasis = getEnrichedBasis(i, -1);
            shp<SparseMatrix> curBasisMatrix(new SparseMatrix(curBasis));
            M_enrichedBasesMatrices.push_back(curBasisMatrix);
            M_enrichedBasesMatricesTransposed.push_back(spcast<SparseMatrix>(curBasisMatrix->transpose()));
        }

        print();
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
        if (M_data("rb/online/basis/useprimalsupremizers", 1))
            msg +=  "\tUsing primal supremizers\n";
        else
            msg +=  "\tNOT using primal supremizers\n";

        if (M_data("rb/online/basis/usedualsupremizers", 0))
        {
            msg +=  "\tUsing dual supremizers\n";
            int nDualSupremizers = M_data("rb/online/basis/numberdualsupremizers", -1);
            if (nDualSupremizers != -1)
                msg +=  "\tUsing " + std::to_string(nDualSupremizers) + " dual supremizers\n";
        }
        else
            msg +=  "\tNOT using dual supremizers\n";

        msg +=  "\tBasis size + supremizers = " + std::to_string(getEnrichedBasis(i,-1).size()) + "\n";
    }
    printlog(WHITE, msg, M_data.getVerbose());
}

void
RBBases::
dump()
{
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
                    std::function<void(shp<VECTOREPETRA>)> transform)
{
    std::vector<shp<VECTOREPETRA>> newEnrichedBasis;
    auto refEnrichedBasis = getEnrichedBasis(index,-1);

    for (auto vec : refEnrichedBasis)
    {
        shp<VECTOREPETRA> newVector(new VECTOREPETRA(M_fespaces[index]->map()));
        *newVector = *vec;

        transform(newVector);
        newEnrichedBasis.push_back(newVector);
    }

    M_enrichedBasesMap[index][ID] = newEnrichedBasis;

    shp<SparseMatrix> basisMatrix(new SparseMatrix(newEnrichedBasis));
    M_enrichedBasesMatricesMap[index][ID] = basisMatrix;
    M_enrichedBasesMatricesTransposedMap[index][ID] = spcast<SparseMatrix>(basisMatrix->transpose());
}

shp<BlockMatrix>
RBBases::
leftProject(shp<BlockMatrix> matrix, unsigned int ID)
{
    shp<BlockMatrix> projectedMatrix(new BlockMatrix(matrix->nRows(), matrix->nCols()));

    for (unsigned int i = 0; i < matrix->nRows(); i++)
    {
        for (unsigned int j = 0; j < matrix->nCols(); j++)
        {
            projectedMatrix->setBlock(i,j,leftProject(spcast<SparseMatrix>(matrix->block(i,j)), i, ID));
        }
    }
    return projectedMatrix;
}

shp<RBMATRIX>
RBBases::
leftProject(shp<FEMATRIX> matrix, unsigned int basisIndex, unsigned int ID)
{
    shp<DenseMatrix> retMat(new DenseMatrix);

    if (matrix->data())
    {
        // the basis is transposed

        auto rangeMap = spcast<MATRIXEPETRA>(getEnrichedBasisMatrices(basisIndex, ID, false)->data())->domainMapPtr();
        auto domainMap = spcast<MATRIXEPETRA>(matrix->data())->domainMapPtr();

        shp<MATRIXEPETRA> innerMatrix(new MATRIXEPETRA(*rangeMap));

        // for some reason it is more efficient to transpose and then multiply
        spcast<MATRIXEPETRA>(getEnrichedBasisMatrices(basisIndex, ID, true)->data())->multiply(false, *spcast<MATRIXEPETRA>(matrix->data()), false,
                                                                        *innerMatrix);

        innerMatrix->globalAssemble(domainMap, rangeMap);

        shp<SparseMatrix> retMatEp(new SparseMatrix());
        retMatEp->setMatrix(innerMatrix);

        return retMatEp->toDenseMatrixPtr();
    }
    return retMat;
}

// we assume that the basis for the lagrange multiplier is the identity
shp<BlockVector>
RBBases::
projectOnLagrangeSpace(shp<BlockVector> vector)
{
    if (vector->nRows() > 1)
        throw new Exception("projectOnLagrangeSpace: error!");

    shp<BlockVector> retVec(new BlockVector(1));
    retVec->setBlock(0,spcast<DistributedVector>(vector->block(0))->toDenseVectorPtr());

    return retVec;
}

shp<BlockVector>
RBBases::
leftProject(shp<BlockVector> vector, unsigned int ID)
{
    shp<BlockVector> projectedVector(new BlockVector(vector->nRows()));

    for (unsigned int i = 0; i < vector->nRows(); i++)
        projectedVector->setBlock(i,leftProject(spcast<DistributedVector>(vector->block(i)), i, ID));

    return projectedVector;
}

shp<DenseVector>
RBBases::
leftProject(shp<DistributedVector> vector, unsigned int basisIndex, unsigned int ID)
{
    if (vector->data())
    {
        shp<DistributedVector> resVector;
        resVector = spcast<DistributedVector>(getEnrichedBasisMatrices(basisIndex, ID, true)->multiplyByVector(vector));

        return resVector->toDenseVectorPtr();
    }
    return nullptr;
}

shp<aVector>
RBBases::
leftProject(shp<aVector> vector, unsigned int basisIndex, unsigned int ID, shp<aMatrix> normMatrix)
{
    if (vector->data())
    {
        auto rangeMap = spcast<MATRIXEPETRA>(getEnrichedBasisMatrices(basisIndex, ID, false)->data())->domainMapPtr();

        shp<VECTOREPETRA> innerVector(new VECTOREPETRA(*rangeMap));
        shp<MATRIXEPETRA> projectionMatrix(new MATRIXEPETRA(*rangeMap));

        auto enrichedBasisIndexEpetra = spcast<MATRIXEPETRA>(getEnrichedBasisMatrices(basisIndex, ID, true)->data());
        if (normMatrix)
        {
            enrichedBasisIndexEpetra->multiply(false, *spcast<MATRIXEPETRA>(normMatrix->data()), false, *projectionMatrix);
            projectionMatrix->multiply(false, *spcast<VECTOREPETRA>(vector->data()), *innerVector);
        }
        else
            enrichedBasisIndexEpetra->multiply(false, *spcast<VECTOREPETRA>(vector->data()), *innerVector);

        innerVector->globalAssemble();

        shp<DistributedVector> retVecEp(new DistributedVector());
        retVecEp->setVector(innerVector);

        return retVecEp->toDenseVectorPtr();
    }
    return nullptr;
}

shp<BlockMatrix>
RBBases::
rightProject(shp<BlockMatrix> matrix, unsigned int ID)
{
    shp<BlockMatrix> projectedMatrix(new BlockMatrix(matrix->nRows(), matrix->nCols()));

    for (unsigned int i = 0; i < matrix->nRows(); i++)
    {
        for (unsigned int j = 0; j < matrix->nCols(); j++)
        {
            projectedMatrix->setBlock(i,j,rightProject(spcast<SparseMatrix>(matrix->block(i,j)), j, ID));
        }
    }
    return projectedMatrix;
}

shp<DenseMatrix>
RBBases::
rightProject(shp<FEMATRIX> matrix, unsigned int basisIndex, unsigned int ID)
{
    shp<DenseMatrix> retMat(new DenseMatrix());

    if (matrix->data())
    {
        // the basis is transposed
        auto rangeMap = spcast<MATRIXEPETRA>(matrix->data())->rangeMapPtr(); ;
        auto domainMap = spcast<MATRIXEPETRA>(getEnrichedBasisMatrices(basisIndex, ID, false)->data())->domainMapPtr();

        shp<MATRIXEPETRA> innerMatrix(new MATRIXEPETRA(*rangeMap));

        spcast<MATRIXEPETRA>(matrix->data())->multiply(false, *spcast<MATRIXEPETRA>(getEnrichedBasisMatrices(basisIndex, ID, false)->data()),
                                false, *innerMatrix);

        innerMatrix->globalAssemble(domainMap, rangeMap);

        shp<SparseMatrix> retMatEp(new SparseMatrix());
        retMatEp->setMatrix(innerMatrix);

        return retMatEp->toDenseMatrixPtr();
    }
    return retMat;
}


std::vector<shp<VECTOREPETRA>>
RBBases::
getBasis(const unsigned int& index)
{
    if (M_onlineTol < 1e-15)
        return M_bases[index];

    std::vector<shp<VECTOREPETRA>> basis;

    for (unsigned int i = 0; i < M_NsOnline[index]; i++)
        basis.push_back(M_bases[index][i]);

    return basis;
}

shp<aMatrix>
RBBases::
matrixProject(shp<aMatrix> matrix,
              unsigned int basisIndexRow,
              unsigned int basisIndexCol,
              unsigned int ID,
              shp<aMatrix> normMatrix)
{
    shp<DenseMatrix> retMat(new DenseMatrix());

    if (matrix->data())
    {
        // Goal: compute Vrow' A Vcol

        // 1. compute A Vcol = AUX
        auto rangeMap = spcast<MATRIXEPETRA>(matrix->data())->rangeMapPtr();
        auto domainMap = spcast<MATRIXEPETRA>(getEnrichedBasisMatrices(basisIndexCol, ID, false)->data())->domainMapPtr();

        shp<MATRIXEPETRA> auxMatrix(new MATRIXEPETRA(*rangeMap));
        auto enrichedBasisColEpetra = spcast<MATRIXEPETRA>(getEnrichedBasisMatrices(basisIndexCol, ID, false)->data());
        spcast<MATRIXEPETRA>(matrix->data())->multiply(false, *enrichedBasisColEpetra,
                                false, *auxMatrix);
        auxMatrix->globalAssemble(domainMap, rangeMap);

        // 2. compute Vrow' AUX
        rangeMap = spcast<MATRIXEPETRA>(getEnrichedBasisMatrices(basisIndexRow, ID, false)->data())->domainMapPtr();

        shp<MATRIXEPETRA> innerMatrix(new MATRIXEPETRA(*rangeMap));
        shp<MATRIXEPETRA> projectionMatrix(new MATRIXEPETRA(*rangeMap));

        auto enrichedBasisIndexEpetra = spcast<MATRIXEPETRA>(getEnrichedBasisMatrices(basisIndexRow, ID, true)->data());
        if (normMatrix)
        {
            enrichedBasisIndexEpetra->multiply(false, *spcast<MATRIXEPETRA>(normMatrix->data()), false, *projectionMatrix);
            projectionMatrix->multiply(false, *auxMatrix, false, *innerMatrix);
        }
        else
            enrichedBasisIndexEpetra->multiply(false, *auxMatrix,false, *innerMatrix);

        innerMatrix->globalAssemble(domainMap, rangeMap);

        shp<SparseMatrix> retMatEp(new SparseMatrix());
        retMatEp->setMatrix(innerMatrix);

        return retMatEp->toDenseMatrixPtr();
    }
    return retMat;
}

std::vector<shp<VECTOREPETRA>>
RBBases::
getEnrichedBasis(const unsigned int& index, unsigned int ID)
{
    if (M_enrichedBasesMap[index][ID].size() > 0)
        return M_enrichedBasesMap[index][ID];

    std::vector<shp<VECTOREPETRA>> retVectors = getBasis(index);
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

shp<FEMATRIX>
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

shp<VECTOREPETRA>
RBBases::
reconstructFEFunction(shp<aVector> rbSolution, unsigned int index,
                      unsigned int ID)
{
    shp<DistributedVector> rbSolutionDistributed;
    if (rbSolution->type() == DENSE)
    {
        rbSolutionDistributed = DistributedVector::convertDenseVector(spcast<DenseVector>(rbSolution),M_comm);
    }
    else
    {
        rbSolutionDistributed = spcast<DistributedVector>(rbSolution);
    }

    shp<DistributedVector> res = spcast<DistributedVector>(getEnrichedBasisMatrices(index, ID, false)->multiplyByVector(rbSolutionDistributed));

    return spcast<VECTOREPETRA>(res->data());
}

// void
// RBBases::
// normalizeBasis(const unsigned int& index, shp<MATRIXEPETRA> normMatrix)
// {
//     throw new Exception("normalizeBasis not implemented in RBBases");
//     // printlog(MAGENTA, "\n[RBBases] normalizing via stabilized Gram Schimdt...\n", M_data.getVerbose());
//     //
//     // const double thrsh = 1e-5;
//     //
//     // auto project = [normMatrix] (const shp<VECTOREPETRA>& vec1,
//     //                              const shp<VECTOREPETRA>& vec2)->double
//     // {
//     //     auto rmap = *normMatrix->rangeMapPtr();
//     //
//     //     VECTOREPETRA aux(rmap);
//     //     normMatrix->matrixPtr()->Multiply(false, vec1->epetraVector(),
//     //                                       aux.epetraVector());
//     //
//     //     return aux.dot(*vec2);
//     // };
//     //
//     // std::vector<bool> keepVector(M_bases[index].size(), true);
//     //
//     // auto orthonormalizeWrtBasis = [&](shp<VECTOREPETRA>& vector,
//     //                                   std::vector<shp<VECTOREPETRA>>& basis,
//     //                                   const unsigned int& count)->void
//     // {
//     //     auto rmap = *normMatrix->rangeMapPtr();
//     //
//     //     double normVector = project(vector, vector);
//     //     *vector /= std::sqrt(normVector);
//     //
//     //     unsigned int basisIndex = 0;
//     //     for (auto& basisV : basis)
//     //     {
//     //         if (keepVector[basisIndex])
//     //         {
//     //             // if this is close to 1 then the vectors are almost parallel
//     //             // because they are both unitary. Note that we assume basisV to be
//     //             // unitary
//     //             double coeff = project(vector, basisV);
//     //
//     //             if (std::abs(1.0 - std::abs(coeff)) > thrsh)
//     //             {
//     //                 shp<VECTOREPETRA> aux(new VECTOREPETRA(rmap));
//     //                 *aux = *vector - (*basisV) * coeff;
//     //
//     //                 vector = aux;
//     //             }
//     //             else
//     //             {
//     //                 // in this case the supremizer is essentially in the column space.
//     //                 // therefore we set it equal to the vector that triggers the check
//     //                 std::string msg;
//     //                 vector = basisV;
//     //                 if (basisIndex > keepVector.size())
//     //                 {
//     //                     msg = "\nAttention: basisIndex > size of keepVector. ";
//     //                     msg += "This indicates that two supremizers are not independent";
//     //                     msg += " between each other\n";
//     //                     printlog(RED, msg, M_data.getVerbose());
//     //                 }
//     //                 else
//     //                 {
//     //                     msg = "\nSwapping supremizer with primal vector\n";
//     //                     printlog(WHITE, msg, M_data.getVerbose());
//     //                     keepVector[basisIndex] = false;
//     //                 }
//     //                 return;
//     //             }
//     //         }
//     //         basisIndex++;
//     //     }
//     //     normVector = project(vector, vector);
//     //     *vector /= std::sqrt(normVector);
//     // };
//     //
//     // // re-orthonormalize the primal basis because if it is not orhonormal to
//     // // machine precision it's a mess
//     //
//     // std::vector<shp<VECTOREPETRA>> incrBasis;
//     //
//     // printlog(GREEN, "normalizing basis\n", M_data.getVerbose());
//     // for (unsigned int i = 0; i < M_bases[index].size(); i++)
//     // {
//     //     // orthonormalizeWrtBasis(M_bases[index][i], incrBasis, 0);
//     //     *M_bases[index][i] /= std::sqrt(project(M_bases[index][i],M_bases[index][i]));
//     //     incrBasis.push_back(M_bases[index][i]);
//     // }
//     //
//     // // double check
//     // // for (unsigned int i = 10; i < 11; i++)
//     // // {
//     // //     for (unsigned int j = 0; j < incrBasis.size(); j++)
//     // //     {
//     // //         double proj = project(incrBasis[i],incrBasis[j]);
//     // //         std::cout << "i = " << i << " j = " << j << " proj = " << proj << std::endl;
//     // //         if (i == j)
//     // //             std::cout << "i = " << i << " j = " << j << " 1 - proj = " << abs(1 - proj) << std::endl;
//     // //     }
//     // // }
//     //
//     // // first orthonormalize primal supremizers
//     // printlog(GREEN, "normalizing primal supremizers\n", M_data.getVerbose());
//     // for (unsigned int j = 0; j < M_numFields; j++)
//     // {
//     //     if (M_primalSupremizers(index,j).size() > 0)
//     //     {
//     //         unsigned int count = 0;
//     //         for (auto& supr : M_primalSupremizers(index,j))
//     //         {
//     //             std::string msg = "orthonormalizing primal supremizer ";
//     //             msg += std::to_string(count);
//     //             msg += "\n";
//     //             printlog(YELLOW, msg, M_data.getVerbose());
//     //
//     //             orthonormalizeWrtBasis(supr, incrBasis, count);
//     //             incrBasis.push_back(supr);
//     //             count++;
//     //         }
//     //     }
//     // }
//     //
//     // // then dual supremizers
//     // printlog(GREEN, "normalizing dual supremizers\n", M_data.getVerbose());
//     // if (M_dualSupremizers[index].size() > 0)
//     // {
//     //     unsigned int count = 0;
//     //     for (auto& supr : M_dualSupremizers[index])
//     //     {
//     //         std::string msg = "orthonormalizing dual supremizer ";
//     //         msg += std::to_string(count);
//     //         msg += "\n";
//     //         printlog(YELLOW, msg, M_data.getVerbose());
//     //
//     //         orthonormalizeWrtBasis(supr, incrBasis, count);
//     //         incrBasis.push_back(supr);
//     //         count++;
//     //     }
//     // }
//     //
//     // // finally select the basis functions that we keep
//     // std::vector<shp<VECTOREPETRA>> newBasis;
//     // for (unsigned int i = 0; i < keepVector.size(); i++)
//     // {
//     //     if (keepVector[i])
//     //         newBasis.push_back(incrBasis[i]);
//     // }
//     //
//     // M_bases[index] = newBasis;
//     //
//     // // // double check
//     // // for (unsigned int i = 0; i < incrBasis.size(); i++)
//     // // {
//     // //     for (unsigned int j = 0; j < incrBasis.size(); j++)
//     // //     {
//     // //         double proj = project(incrBasis[i],incrBasis[j]);
//     // //         std::cout << "i = " << i << " j = " << j << " proj = " << proj << std::endl;
//     // //     }
//     // // }
// }

}  // namespace RedMA
