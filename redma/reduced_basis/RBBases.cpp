#include "RBBases.hpp"

namespace RedMA
{

RBBases::
RBBases(const DataContainer& data, EPETRACOMM comm) :
  M_data(data),
  M_comm(comm)
{

}

void
RBBases::
setNumberOfFields(const unsigned int& nfields)
{
    M_numFields = nfields;
    M_bases.resize(nfields);
    M_svs.resize(nfields);
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
}

void
RBBases::
dump()
{
    using namespace boost::filesystem;

    bool binary = M_data("rbbasis/dumpbinary", true);

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

std::vector<SHP(VECTOREPETRA)>&
RBBases::
getBasis(const unsigned int& index, double tol)
{
    return M_bases[index];
}

std::vector<SHP(VECTOREPETRA)>
RBBases::
getEnrichedBasis(const unsigned int& index, double tol)
{
    std::vector<SHP(VECTOREPETRA)> retVectors = getBasis(index, tol);

    for (unsigned int j = 0; j < M_numFields; j++)
    {
        for (auto primal : M_primalSupremizers(index,j))
            retVectors.push_back(primal);
    }

    for (auto dual : M_dualSupremizers[index])
        retVectors.push_back(dual);

    return retVectors;
}
}  // namespace RedMA
