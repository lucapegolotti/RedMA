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

void
AbstractAssembler::
addDualMaps(MapEpetraPtr& globalMap, std::vector<unsigned int>& dimensions)
{
    for (MapVectorSTD::iterator it = M_dualMaps.begin();
         it != M_dualMaps.end(); it++)
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

void
AbstractAssembler::
buildLagrangeMultiplierBasisFourier(const unsigned int& frequencies,
                                    const unsigned int& nComponents,
                                    ETFESpaceCouplingPtr couplingFespace,
                                    MapEpetraPtr primalMap,
                                    FESpacePtr primalFespace,
                                    GeometricFace face,
                                    const unsigned int& faceFlag)
{
    using namespace LifeV;
    using namespace ExpressionAssembly;

    const Real dropTolerance(2.0 * std::numeric_limits<Real>::min());
    QuadratureBoundary boundaryQuadRule(buildTetraBDQR(quadRuleTria7pt));

    std::shared_ptr<FourierBasisFunction> basisFunction;
    basisFunction.reset(new FourierBasisFunction(face, frequencies));

    unsigned int nLagMultiplierBFs = basisFunction->getNumBasisFunctions();
    VectorPtr* couplingVectors = new VectorPtr[nLagMultiplierBFs];

    MapEpetra couplingMap = couplingFespace->map();
    MatrixPtr boundaryMassMatrix(new Matrix(couplingMap));

    MeshPtr mesh = couplingFespace->mesh();

    // assemble boundary mass matrix to orthonormalize w.r.t. L2 product
    integrate(boundary(mesh, faceFlag),
              boundaryQuadRule,
              couplingFespace,
              couplingFespace,
              phi_i * phi_j
              ) >> boundaryMassMatrix;

    for (unsigned int i = 0; i < nLagMultiplierBFs; i++)
    {
        VectorPtr currentMode(new Vector(couplingMap, Repeated));

        basisFunction->setIndex(i);

        integrate(boundary(mesh, faceFlag),
                  boundaryQuadRule,
                  couplingFespace,
                  eval(basisFunction, X) * phi_i
              ) >> currentMode;

        couplingVectors[i] = currentMode;
    }
    gramSchmidt(couplingVectors, boundaryMassMatrix, nLagMultiplierBFs);

    unsigned int myel = (nComponents * nLagMultiplierBFs) / M_comm->NumProc();

    // the first process takes care of the remainder
    if (M_comm->MyPID() == 0)
    {
        myel += (nComponents * nLagMultiplierBFs) % M_comm->NumProc();
    }

    MapEpetraPtr lagrangeMultiplierMap;
    lagrangeMultiplierMap.reset(new MapEpetra(nComponents * nLagMultiplierBFs,
                                              myel, 0, M_comm));
    M_dualMaps.push_back(lagrangeMultiplierMap);

    // note:: we have to specify the second argument of the constructor (number
    // of elements per row)
    MatrixPtr QT(new Matrix(*primalMap, nLagMultiplierBFs, false));
    QT->zero();

    MatrixPtr Q(new Matrix(*lagrangeMultiplierMap));
    Q->zero();

    Epetra_Map primalMapEpetra = couplingVectors[0]->epetraMap();
    unsigned int numElements = primalMapEpetra.NumMyElements();
    unsigned int nTotalDofs = primalFespace->dof().numTotalDof();
    for (unsigned int dim = 0; dim < nComponents; dim++)
    {
        for (unsigned int i = 0; i < nLagMultiplierBFs; i++)
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
                        QT->addToCoefficient(gdof + dim * nTotalDofs,
                                             i + dim * nLagMultiplierBFs,
                                             value);
                    }
                }
            }
        }
    }
    QT->globalAssemble(lagrangeMultiplierMap, primalMap);

    Q = QT->transpose();
    Q->globalAssemble(primalMap, lagrangeMultiplierMap);

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
gramSchmidt(VectorPtr* basis, MatrixPtr massMatrix,
            unsigned int& nVectors)
{
    for (unsigned int i = 0; i < nVectors; i++)
    {
        VectorPtr v(new Vector(*basis[i], LifeV::Unique));
        basis[i] = v;
    }

    bool* quickSearch = new bool[nVectors];
    unsigned int vectorsZeroed = 0;
    double* norms = new double[nVectors];
    for (unsigned int i = 0; i < nVectors; i++)
    {
        double u1 = std::sqrt(dotProd(basis[i], basis[i], massMatrix));
        norms[i] = u1;
        quickSearch[i] = false;
        if (u1 < 5e-12)
        {
            quickSearch[i] = true;
            basis[i]->zero();
            vectorsZeroed++;
        }
        else
        {
            for (unsigned int j = 0; j < i; j++)
            {
                if (!quickSearch[j])
                {
                    double u2 = dotProd(basis[i], basis[j], massMatrix);
                    double theta = std::acos(std::abs(u2)/(norms[i]*norms[j]));

                    if (std::abs(theta) < 5e-1)
                    {
                        quickSearch[i] = true;
                        basis[i]->zero();
                        vectorsZeroed++;
                        break;
                    }
                }
            }
        }
    }
    unsigned int auxIndex = 0;

    std::string msg("[AbstractAssembler] GramSchmidt -> ");
    msg += "reducing number of bfs from " + std::to_string(nVectors) + " to" +
           std::to_string(nVectors - vectorsZeroed) + "\n";
    printlog(MAGENTA, msg, M_verbose);

    for (unsigned int i = 0; i < nVectors - vectorsZeroed; i++)
    {
        while (quickSearch[auxIndex] && auxIndex < nVectors)
        {
            auxIndex++;
        }
        basis[i] = basis[auxIndex];
    }

    nVectors -= vectorsZeroed;

    // normalize vectors. WARNING: it might be that this must be done on the
    // whole vector (comprising also the adjacent domain)
    for (unsigned int i = 0; i < nVectors; i++)
    {
        for (unsigned int j = 0; j < i; j++)
        {
            double uv = dotProd(basis[i], basis[j], massMatrix);
            double uu = dotProd(basis[j], basis[j], massMatrix);

            *basis[i] += (-uv/uu) * (*basis[j]);
        }

        double curnorm = dotProd(basis[i], basis[i], massMatrix);
        *basis[i] *= (1.0/std::sqrt(curnorm));
    }

    delete[] quickSearch;
    delete[] norms;
}

double
AbstractAssembler::
dotProd(VectorPtr vector1, VectorPtr vector2, MatrixPtr mass)
{
    VectorPtr aux(new Vector(*vector1));
    *aux = (*mass) * (*vector1);
    double res = 0;
    aux->dot(*vector2, res);
    return res;
}

}  // namespace RedMA
