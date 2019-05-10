#include <AbstractAssembler.hpp>

namespace RedMA
{

AbstractAssembler::
AbstractAssembler(const GetPot& datafile, commPtr_Type comm,
                  const TreeNodePtr& treeNode, bool verbose) :
  M_datafile(datafile),
  M_comm(comm),
  M_treeNode(treeNode),
  M_verbose(verbose)
{
}

void
AbstractAssembler::
addPrimalMaps(MapEpetraPtr& globalMap, std::vector<unsigned int>& dimensions)
{
    for (MapVectorSTD::iterator it = M_primalMaps.begin();
         it != M_primalMaps.end(); it++)
    {
        *globalMap += *(*it);
        dimensions.push_back((*it)->map(LifeV::Unique)->NumGlobalElements());
    }
}

std::vector<AbstractAssembler::MapEpetraPtr>
AbstractAssembler::
getPrimalMapVector()
{
    return M_primalMaps;
}

std::vector<AbstractAssembler::MapEpetraPtr>
AbstractAssembler::
getDualMapVector()
{
    return M_dualMaps;
}

void
AbstractAssembler::
gramSchmidt(AbstractAssembler::VectorPtr* basis1,
            AbstractAssembler::MatrixPtr massMatrix1,
            AbstractAssembler::VectorPtr* basis2,
            AbstractAssembler::MatrixPtr massMatrix2,
            unsigned int& nVectors)
{
    for (unsigned int i = 0; i < nVectors; i++)
    {
        VectorPtr v1(new Vector(*basis1[i], LifeV::Unique));
        VectorPtr v2(new Vector(*basis2[i], LifeV::Unique));
        basis1[i] = v1;
        basis2[i] = v2;
    }

    bool* quickSearch = new bool[nVectors];
    unsigned int vectorsZeroed = 0;
    double* norms = new double[nVectors];
    for (unsigned int i = 0; i < nVectors; i++)
    {
        double u1 = std::sqrt(dotProd(basis1, basis2, i, i,
                                      massMatrix1, massMatrix2));
        norms[i] = u1;
        quickSearch[i] = false;
        if (u1 < 5e-12)
        {
            quickSearch[i] = true;
            basis1[i]->zero();
            basis2[i]->zero();
            vectorsZeroed++;
        }
        else
        {
            for (unsigned int j = 0; j < i; j++)
            {
                if (!quickSearch[j])
                {
                    double u2 = dotProd(basis1, basis2, i, j,
                                        massMatrix1, massMatrix2);
                    double theta = std::acos(std::abs(u2)/(norms[i]*norms[j]));

                    if (std::abs(theta) < 5e-1)
                    {
                        quickSearch[i] = true;
                        basis1[i]->zero();
                        basis2[i]->zero();
                        vectorsZeroed++;
                        break;
                    }
                }
            }
        }
    }
    unsigned int auxIndex = 0;

    for (unsigned int i = 0; i < nVectors - vectorsZeroed; i++)
    {
        while (quickSearch[auxIndex] && auxIndex < nVectors)
        {
            auxIndex++;
        }
        basis1[i] = basis1[auxIndex];
        basis2[i] = basis2[auxIndex];
        auxIndex++;
    }

    nVectors -= vectorsZeroed;

    // normalize vectors. WARNING: it might be that this must be done on the
    // whole vector (comprising also the adjacent domain)
    for (unsigned int i = 0; i < nVectors; i++)
    {
        for (unsigned int j = 0; j < i; j++)
        {
            double uv = AbstractAssembler::dotProd(basis1, basis2, i, j,
                                                   massMatrix1, massMatrix2);
            double uu = AbstractAssembler::dotProd(basis1, basis2, j, j,
                                                   massMatrix1, massMatrix2);

            *basis1[i] += (-uv/uu) * (*basis1[j]);
            *basis2[i] += (-uv/uu) * (*basis2[j]);
        }

        double curnorm = AbstractAssembler::dotProd(basis1, basis2, i, i,
                                                    massMatrix1, massMatrix2);
        *basis1[i] *= (1.0/std::sqrt(curnorm));
        *basis2[i] *= (1.0/std::sqrt(curnorm));
    }

    delete[] quickSearch;
    delete[] norms;
}

double
AbstractAssembler::
dotProd(VectorPtr* basis1, VectorPtr* basis2, unsigned int index1,
        unsigned int index2, MatrixPtr mass1, MatrixPtr mass2)
{
    double prod1 = 0, prod2 = 0;
    VectorPtr aux1(new Vector(basis1[index1]->map()));
    VectorPtr aux2(new Vector(basis2[index2]->map()));
    *aux1 = (*mass1) * (*basis1[index1]);
    *aux2 = (*mass2) * (*basis2[index2]);
    aux1->dot(*basis1[index1], prod1);
    aux2->dot(*basis2[index2], prod2);
    return prod1 + prod2;
}

AbstractAssembler::VectorPtr*
AbstractAssembler::
assembleCouplingVectorsFourier(const unsigned int& frequencies,
                               const unsigned int& nBasisFunctions,
                               GeometricFace face, const double& coeff)
{
    using namespace LifeV;
    using namespace ExpressionAssembly;

    QuadratureBoundary boundaryQuadRule(buildTetraBDQR(quadRuleTria7pt));

    std::shared_ptr<FourierBasisFunction> basisFunction;
    basisFunction.reset(new FourierBasisFunction(face, frequencies));

    VectorPtr* couplingVectors = new VectorPtr[nBasisFunctions];
    MapEpetra couplingMap = M_couplingFESpaceETA->map();
    MatrixPtr boundaryMassMatrix(new Matrix(couplingMap));

    unsigned int faceFlag = face.M_flag;
    MeshPtr mesh = M_couplingFESpaceETA->mesh();

    for (unsigned int i = 0; i < nBasisFunctions; i++)
    {
        VectorPtr currentMode(new Vector(couplingMap, Repeated));

        basisFunction->setIndex(i);
        integrate(boundary(mesh, faceFlag),
                  boundaryQuadRule,
                  M_couplingFESpaceETA,
                  value(coeff) * eval(basisFunction, X) * phi_i
              ) >> currentMode;
        couplingVectors[i] = currentMode;
    }
    return couplingVectors;
}

AbstractAssembler::MatrixPtr
AbstractAssembler::
assembleBoundaryMatrix(GeometricFace face)
{
    using namespace LifeV;
    using namespace ExpressionAssembly;

    QuadratureBoundary boundaryQuadRule(buildTetraBDQR(quadRuleTria7pt));

    unsigned int faceFlag = face.M_flag;
    MeshPtr mesh = M_couplingFESpaceETA->mesh();

    MapEpetra couplingMap = M_couplingFESpaceETA->map();
    MatrixPtr boundaryMassMatrix(new Matrix(couplingMap));

    // assemble boundary mass matrix to orthonormalize w.r.t. L2 product
    integrate(boundary(mesh, faceFlag),
              boundaryQuadRule,
              M_couplingFESpaceETA,
              M_couplingFESpaceETA,
              phi_i * phi_j
              ) >> boundaryMassMatrix;
    boundaryMassMatrix->globalAssemble();

    return boundaryMassMatrix;
}

void
AbstractAssembler::
fillMatricesWithVectors(VectorPtr* couplingVectors,
                        const unsigned int& nBasisFunctions,
                        MapEpetraPtr lagrangeMap,
                        const unsigned int& flagAdjacentDomain)
{
    using namespace LifeV;
    const Real dropTolerance(2.0 * std::numeric_limits<Real>::min());
    unsigned faceFlag = flagAdjacentDomain;
    // note:: we have to specify the second argument of the constructor (number
    // of elements per row)
    MatrixPtr QT(new Matrix(*M_primalMaps[M_indexCoupling],
                            nBasisFunctions, false));
    QT->zero();

    MatrixPtr Q(new Matrix(*lagrangeMap));
    Q->zero();

    Epetra_Map primalMapEpetra = couplingVectors[0]->epetraMap();
    unsigned int numElements = primalMapEpetra.NumMyElements();
    // unsigned int nTotalDofs = primalFespace->dof().numTotalDof();
    unsigned int nTotalDofs = couplingVectors[0]->size();
    for (unsigned int dim = 0; dim < numberOfComponents(); dim++)
    {
        for (unsigned int i = 0; i < nBasisFunctions; i++)
        {
            Vector couplingVectorUnique(*couplingVectors[i], Unique);
            for (unsigned int dof = 0; dof < numElements; dof++)
            {
                unsigned int gdof = primalMapEpetra.GID(dof);
                if (couplingVectorUnique.isGlobalIDPresent(gdof))
                {
                    double value(couplingVectorUnique[gdof]);
                    if (std::abs(value) > dropTolerance)
                    {
                        // this must be changed to insert coefficients
                        QT->addToCoefficient(gdof + dim * nTotalDofs,
                                             i + dim * nBasisFunctions,
                                             value);
                    }
                }
            }
        }
    }
    QT->globalAssemble(lagrangeMap, M_primalMaps[M_indexCoupling]);

    Q = QT->transpose();
    Q->globalAssemble(M_primalMaps[M_indexCoupling], lagrangeMap);

    if (M_mapQTs.find(faceFlag) == M_mapQTs.end() &&
        M_mapQs.find(faceFlag) == M_mapQs.end())
    {
        M_mapQTs[faceFlag] = QT;
        M_mapQs[faceFlag] = Q;
    }
    else
    {
        std::string errorMsg("Coupling matrices with key = ");
        errorMsg += std::to_string(faceFlag) + " have already been assembled!\n";

        throw Exception(errorMsg);
    }

    delete[] couplingVectors;
}

void
AbstractAssembler::
assembleCouplingMatrices(AbstractAssembler& child,
                         const unsigned int& indexOutlet,
                         const unsigned int& interfaceIndex,
                         AbstractAssembler::MapEpetraPtr& globalMap,
                         std::vector<unsigned int>& dimensions)
{
    std::string msg("[AbstractAssembler] start ");
    msg += "building coupling matrices between blocks with indices ";
    msg += std::to_string(M_treeNode->M_ID);
    msg += " and ";
    msg += std::to_string(child.M_treeNode->M_ID);
    msg += "\n";
    printlog(MAGENTA, msg, M_verbose);

    unsigned int nComponents = numberOfComponents();
    unsigned int nBasisFunctions;
    MapEpetraPtr lagrangeMultiplierMap;

    std::string typeBasis = M_datafile("coupling/type", "fourier");

    if (!std::strcmp(typeBasis.c_str(), "fourier"))
    {
        unsigned int frequencies = M_datafile("coupling/frequencies", 1);
        nBasisFunctions = (2 * frequencies + 1) * (frequencies + 1);

        GeometricFace outlet = M_treeNode->M_block->getOutlet(indexOutlet);
        VectorPtr* couplingVectorsFather =
                    assembleCouplingVectorsFourier(frequencies, nBasisFunctions,
                                                   outlet, 1);
        MatrixPtr massMatrixFather = assembleBoundaryMatrix(outlet);

        GeometricFace inlet = child.M_treeNode->M_block->getInlet();
        VectorPtr* couplingVectorsChild =
              child.assembleCouplingVectorsFourier(frequencies, nBasisFunctions,
                                                     inlet, -1);

        MatrixPtr massMatrixChild = child.assembleBoundaryMatrix(inlet);
        unsigned int prev = nBasisFunctions;
        gramSchmidt(couplingVectorsFather, massMatrixFather, couplingVectorsChild,
                    massMatrixChild, nBasisFunctions);

        msg = "GramSchmidt -> ";
        msg += "reducing number of bfs from " + std::to_string(prev) + " to " +
               std::to_string(nBasisFunctions) + "\n";
        printlog(GREEN, msg, M_verbose);

        // build map for the lagrange multipliers (nBasisFunctions has been
        // modified in gramSchmidt)
        unsigned int myel = (nComponents * nBasisFunctions) / M_comm->NumProc();
        // the first process takes care of the remainder
        if (M_comm->MyPID() == 0)
        {
            myel += (nComponents * nBasisFunctions) % M_comm->NumProc();
        }

        unsigned int mapSize = nComponents * nBasisFunctions;
        lagrangeMultiplierMap.reset(new
                        AbstractAssembler::MapEpetra(mapSize, myel,
                                                     0, M_comm));
        M_dualMaps.push_back(lagrangeMultiplierMap);
        child.M_dualMaps.push_back(lagrangeMultiplierMap);

        fillMatricesWithVectors(couplingVectorsFather, nBasisFunctions,
                                lagrangeMultiplierMap,
                                child.M_treeNode->M_ID);

        child.fillMatricesWithVectors(couplingVectorsChild, nBasisFunctions,
                                      lagrangeMultiplierMap, M_treeNode->M_ID);
    }

    *globalMap += *lagrangeMultiplierMap;
    dimensions.push_back(lagrangeMultiplierMap->map(LifeV::Unique)
                                              ->NumGlobalElements());
    M_interfacesIndices.push_back(interfaceIndex);
    child.M_interfacesIndices.push_back(interfaceIndex);
    printlog(MAGENTA, "done\n", M_verbose);
}

// we copy the matrix in order to control boundary conditions
AbstractAssembler::MatrixPtr
AbstractAssembler::
getQT(const unsigned int& flag)
{
    MatrixPtr retMatrix(new Matrix(*M_mapQTs[flag]));
    return retMatrix;
}

AbstractAssembler::MatrixPtr
AbstractAssembler::
getQ(const unsigned int& flag)
{
    MatrixPtr retMatrix(new Matrix(*M_mapQs[flag]));
    return retMatrix;
}

unsigned int
AbstractAssembler::
getIndexCoupling()
{
    return M_indexCoupling;
}

std::vector<unsigned int>
AbstractAssembler::
getInterfacesIndices()
{
    return M_interfacesIndices;
}

}  // namespace RedMA
