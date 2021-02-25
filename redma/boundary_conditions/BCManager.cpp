#include "BCManager.hpp"

namespace RedMA
{

BCManager::
BCManager(const DataContainer& data, shp<TreeNode> treeNode) :
  M_data(data),
  M_treeNode(treeNode)
{
    M_inletBCType = data("bc_conditions/inlet_bc_type", "dirichlet");

    M_inflow = data.getInflow();
    M_strongDirichlet = std::strcmp(data("bc_conditions/inletdirichlet", "weak").c_str(),"strong") == 0;
    M_coefficientInflow = data("bc_conditions/coefficientinflow", 1.0);

    if (M_treeNode->isOutletNode())
        this->parseOutflowNeumannData();
    if (M_treeNode->isInletNode())
        this->checkInflowLaw();

    M_inletFlag = treeNode->M_block->getInlet().M_flag;
    for(auto out_face : treeNode->M_block->getOutlets())
        M_outletFlags.push_back(out_face.M_flag);
    M_wallFlag = treeNode->M_block->getWallFlag();
    M_inletRingFlag = treeNode->M_block->getInlet().M_ringFlag;
    M_outletRingFlag = treeNode->M_block->getOutlet(0).M_ringFlag;  // assuming all outlets are the same
}

void
BCManager::
parseOutflowNeumannData()
{
    unsigned int numConditions = M_data("bc_conditions/numoutletbcs", 0);

    for (unsigned int outletIndex = 0; outletIndex < numConditions; outletIndex++)
    {
        std::string dataEntry = "bc_conditions/outlet" + std::to_string(outletIndex);

        unsigned int blockindex = M_data(dataEntry + "/blockindex", 0);
        if (M_treeNode->M_ID == blockindex)
        {
            unsigned int boundaryflag = M_data(dataEntry + "/boundaryflag", 2);
            M_models[boundaryflag].reset(new WindkesselModel(M_data, dataEntry, outletIndex));
        }
    }
}

void
BCManager::
applyInflowDirichletBCs(shp<LifeV::BCHandler> bcs, const bool& zeroFlag) const
{
    std::function<double(double)> inflowLaw;
    if (zeroFlag)
        inflowLaw = std::function<double(double)>([](double t){return 0.0;});
    else
        inflowLaw = this->M_inflow;

    auto foo = std::bind(poiseuilleInflow,
                         std::placeholders::_1,
                         std::placeholders::_2,
                         std::placeholders::_3,
                         std::placeholders::_4,
                         std::placeholders::_5,
                         M_treeNode->M_block->getInlet(),
                         inflowLaw,
                         this->M_coefficientInflow);

    LifeV::BCFunctionBase inflowFunction(foo);
    bcs->addBC("Inlet", M_inletFlag, LifeV::Essential, LifeV::Full,
               inflowFunction, 3);
}

void
BCManager::
applyInflowNeumannBCs(shp<LifeV::BCHandler> bcs, const bool& zeroFlag) const
{
    std::function<double(double)> inflowLaw;
    if (zeroFlag)
        inflowLaw = std::function<double(double)>([](double t){return 0.0;});
    else
        inflowLaw = this->M_inflow;

    auto inflowBoundaryCondition = std::bind(neumannInflow,
                                             std::placeholders::_1,
                                             std::placeholders::_2,
                                             std::placeholders::_3,
                                             std::placeholders::_4,
                                             std::placeholders::_5,
                                             inflowLaw);

    LifeV::BCFunctionBase inflowFunction(inflowBoundaryCondition);
    bcs->addBC("Inlet", M_inletFlag, LifeV::Natural, LifeV::Normal,
               inflowFunction);
}

void
BCManager::
addInletBC(shp<LifeV::BCHandler> bcs, const bool& ringOnly, const bool& zeroFlag) const
{
    if (M_treeNode->isInletNode())
    {
        LifeV::BCFunctionBase zeroFunction(fZero);

        if (!std::strcmp(M_inletBCType.c_str(), "dirichlet"))
            // imposing an inflow velocity profile
            this->applyInflowDirichletBCs(bcs, zeroFlag);

        else if (!std::strcmp(M_inletBCType.c_str(), "neumann"))
            // imposing an inflow pressure --> actually this is done at RHS so this is not called!
            this->applyInflowNeumannBCs(bcs, zeroFlag);

        else
            throw new Exception("Selected inflow BCs not implemented!");

        /*if (!ringOnly) {
            bcs->addBC("Wall", M_wallFlag, LifeV::Essential,
                       LifeV::Full, zeroFunction, 3);
        }
        else {
            bcs->addBC("InletRing", M_inletRingFlag, LifeV::EssentialEdges,
                       LifeV::Full, zeroFunction, 3);

            *//*std::vector <LifeV::ID> compz (1);
            compz[0] = 2;
            bcs->addBC("InletRingZ", M_inletRingFlag, LifeV::EssentialEdges,
                       LifeV::Component, zeroFunction, compz);*//*
        }*/
    }
}

void
BCManager::
applyDirichletBCs(const double& time, BlockVector& input,
                  shp<FESPACE> fespace, const unsigned int& index,
                  const bool& ringOnly)
{
    shp<LifeV::BCHandler> bcs = createBCHandler0Dirichlet(ringOnly);

    if (!std::strcmp(M_inletBCType.c_str(), "dirichlet"))
        addInletBC(bcs, ringOnly);

    bcs->bcUpdate(*fespace->mesh(), fespace->feBd(), fespace->dof());

    shp<VECTOREPETRA> curVec(spcast<VECTOREPETRA>(input.block(index)->data()));

    if (curVec)
        bcManageRhs(*curVec, *fespace->mesh(), fespace->dof(),
                    *bcs, fespace->feBd(), 1.0, time);

    if (ringOnly)
    {
        std::cout<<"DBC: 1"<<std::endl;

        shp<LifeV::BCHandler> bcsRing = createBCHandler0DirichletRing();
        bcsRing->bcUpdate(*fespace->mesh(), fespace->feBd(), fespace->dof());

        std::cout<<"DBC: 2"<<std::endl;

        // shift to (t1,t2,n) reference system
        this->shiftToNormalTangentialCoordSystem(nullptr, curVec, fespace);

        std::cout<<"DBC: 3"<<std::endl;

        // apply proper ring BCs
        bcManageRhs(*curVec, *fespace->mesh(), fespace->dof(),
                *bcsRing, fespace->feBd(), 0.0, 0.0);

        std::cout<<"DBC: 4"<<std::endl;

        // shift back to (x,y,z) reference system
        this->shiftToCartesianCoordSystem(nullptr, curVec, fespace);

        std::cout<<"DBC: 5"<<std::endl;
    }
}

void
BCManager::
apply0DirichletMatrix(BlockMatrix& input,
                      shp<FESPACE> fespace,
                      const unsigned int& index,
                      const double& diagCoefficient,
                      const bool& ringOnly)
{
    shp<LifeV::BCHandler> bcs = createBCHandler0Dirichlet(ringOnly);

    if ((!std::strcmp(M_inletBCType.c_str(), "dirichlet")) && (M_strongDirichlet))
        addInletBC(bcs, ringOnly, true);

    bcs->bcUpdate(*fespace->mesh(), fespace->feBd(), fespace->dof());

    unsigned int nCols = input.nCols();
    for (unsigned int j = 0; j < nCols; j++)
    {
        if (!input.block(index, j)->isZero())
        {
            shp<MATRIXEPETRA> curMatrix = spcast<MATRIXEPETRA>(input.block(index, j)->data());
            auto domainMap = curMatrix->domainMapPtr();
            auto rangeMap = curMatrix->rangeMapPtr();

            bcManageMatrix(*curMatrix, *fespace->mesh(),
                           fespace->dof(), *bcs, fespace->feBd(),
                           (j == index) * diagCoefficient, 0.0);

            if (ringOnly && !index && !j)
            {
                std::cout<<"DBC0_mat: 1"<<std::endl;

                // shift to (t1,t2,n) reference system
                this->shiftToNormalTangentialCoordSystem(curMatrix, nullptr, fespace);

                std::cout<<"DBC0_mat: 2"<<std::endl;

                // apply BCs
                bcManageMatrix(*curMatrix, *fespace->mesh(),
                               fespace->dof(), *bcs, fespace->feBd(),
                               (j == index) * diagCoefficient, 0.0);

                std::cout<<"DBC0_mat: 3"<<std::endl;

                // shift back to in (x,y,z) reference system
                this->shiftToCartesianCoordSystem(curMatrix, nullptr, fespace);

                std::cout<<"DBC0_mat: 4"<<std::endl;
            }

            curMatrix->globalAssemble(domainMap, rangeMap);
        }
    }
}

void
BCManager::
apply0DirichletBCs(BlockVector& input, shp<FESPACE> fespace,
                   const unsigned int& index,
                   const bool& ringOnly)
{
    shp<LifeV::BCHandler> bcs = createBCHandler0Dirichlet(ringOnly);

    if ((!std::strcmp(M_inletBCType.c_str(), "dirichlet")) && (M_strongDirichlet))
        addInletBC(bcs, ringOnly, true);

    bcs->bcUpdate(*fespace->mesh(), fespace->feBd(), fespace->dof());

    shp<VECTOREPETRA> curVec(spcast<VECTOREPETRA>(input.block(index)->data()));
    if (curVec)
        bcManageRhs(*curVec, *fespace->mesh(), fespace->dof(),
                    *bcs, fespace->feBd(), 0.0, 0.0);

    if (ringOnly)
    {
        std::cout<<"DBC0: 1"<<std::endl;

        shp<LifeV::BCHandler> bcsRing = createBCHandler0DirichletRing();
        bcsRing->bcUpdate(*fespace->mesh(), fespace->feBd(), fespace->dof());

        std::cout<<"DBC0: 2"<<std::endl;

        // shift to (t1,t2,n) reference system
        this->shiftToNormalTangentialCoordSystem(nullptr, curVec, fespace);

        std::cout<<"DBC0: 3"<<std::endl;

        // apply proper ring BCs
        bcManageRhs(*curVec, *fespace->mesh(), fespace->dof(),
                    *bcsRing, fespace->feBd(), 0.0, 0.0);

        std::cout<<"DBC0: 4"<<std::endl;

        // shift back to (x,y,z) reference system
        this->shiftToCartesianCoordSystem(nullptr, curVec, fespace);

        std::cout<<"DBC0: 5"<<std::endl;
    }
}

shp<LifeV::BCHandler>
BCManager::
createBCHandler0Dirichlet(const bool& ringOnly) const
{
    LifeV::BCFunctionBase zeroFunction(fZero);

    shp<LifeV::BCHandler> bcs;
    bcs.reset(new LifeV::BCHandler);

    if (!ringOnly)
        bcs->addBC("Wall", M_wallFlag, LifeV::Essential,
                   LifeV::Full, zeroFunction, 3);
    /*else {
        *//*std::vector <LifeV::ID> compz (1);
        compz[0] = 2;*//*

        if (M_treeNode->isInletNode()) {
             bcs->addBC("InletRing", M_inletRingFlag, LifeV::EssentialEdges,
                          LifeV::Full, zeroFunction, 3);

            // TODO: generalize for a generic inlet orientation
            *//*bcs->addBC("InletRingZ", M_inletRingFlag, LifeV::EssentialEdges,
                       LifeV::Component, zeroFunction, compz);*//*
        }

        if (M_treeNode->isOutletNode()) {
            // this imposes the BC at all outlets, as the outlet ring flag is the same
            // for all outlets; thus outlet building blocks MUST be tubes!!

            bcs->addBC("OutletRing", M_outletRingFlag, LifeV::EssentialEdges,
                         LifeV::Full, zeroFunction, 3);

            // TODO: generalize for a generic outlet orientation
            *//*bcs->addBC("OutletRingZ", M_outletRingFlag, LifeV::EssentialEdges,
                       LifeV::Component, zeroFunction, compz);*//*
        }
    }*/

    return bcs;
}

shp<LifeV::BCHandler>
BCManager::
createBCHandler0DirichletRing() const
{
    LifeV::BCFunctionBase zeroFunction(fZero);

    shp<LifeV::BCHandler> bcs;
    bcs.reset(new LifeV::BCHandler);

    std::vector <LifeV::ID> compz (1);
        compz[0] = 2;

    if (M_treeNode->isInletNode()) {
        /*bcs->addBC("InletRing", M_inletRingFlag, LifeV::EssentialEdges,
                   LifeV::Full, zeroFunction, 3);*/

        bcs->addBC("InletRingZ", M_inletRingFlag, LifeV::EssentialEdges,
                   LifeV::Component, zeroFunction, compz);
    }

    if (M_treeNode->isOutletNode()) {
        /*bcs->addBC("OutletRing", M_outletRingFlag, LifeV::EssentialEdges,
                   LifeV::Full, zeroFunction, 3);*/

        bcs->addBC("OutletRingZ", M_outletRingFlag, LifeV::EssentialEdges,
                   LifeV::Component, zeroFunction, compz);
    }

    return bcs;
}

double
BCManager::
poiseuilleInflow(const double& t, const double& x, const double& y,
                 const double& z, const unsigned int& i,
                 const GeometricFace& face, const std::function<double(double)> inflow,
                 const double& coefficient)
{
    const Vector3D& center = face.M_center;
    const Vector3D& normal = face.M_normal;
    double R = face.M_radius;

    Vector3D curPoint(x,y,z);
    Vector3D diff = curPoint - center;
    double r = diff.norm();

    // we suppose that inflow is the flowrate and we want to find the max velocity
    const double maxU = inflow(t) * 2.0 / (M_PI * R * R);
    double inflowNorm = maxU * (1.0 - (r * r)/(R * R));

    if (inflow(t) < 0)
        inflowNorm = inflowNorm < 0 ? inflowNorm : 0;

    inflowNorm = inflowNorm > 0 ? inflowNorm : 0;

    Vector3D inflowValue = -inflowNorm * normal * coefficient;
    return inflowValue[i];
}

double
BCManager::
neumannInflow(const double &t, const double &x,
              const double &y, const double &z, const unsigned int &i,
              std::function<double(double)> inflowLaw)
{
    return inflowLaw(t);
}

void
BCManager::
checkInflowLaw()
{
    if ( (!M_inflow) and (std::strcmp(M_inletBCType.c_str(), "dirichlet") == 0) ) {
        printlog(RED, "Setting default Dirichlet inflow BC (unit constant))",
                 this->M_data.getVerbose());

        std::function<double(double)> inflow(
                [](double t) {
                    return 1.0;
                });

        M_inflow = inflow;
    }

    else if ( (!M_inflow) and (std::strcmp(M_inletBCType.c_str(), "neumann") == 0) ) {
        printlog(RED, "Setting default Neumann inflow BC (pressure bump at initial times)",
                 this->M_data.getVerbose());

        const double T = 3e-3;
        const double omega = 2.0 * M_PI / T;
        const double Pmax = 13300.0;

        std::function<double(double)> inflow(
                [T, omega, Pmax](double t) {
                    if (t <= T) {
                        return -0.5 * (1.0 - std::cos(omega * t)) * Pmax;
                    }
                    return 0.0;
                });

        M_inflow = inflow;
    }

}

double
BCManager::
getOutflowNeumannBC(const double& time, const double& flag, const double& rate)
{
    auto it = M_models.find(flag);
    if (it == M_models.end())
        return 0.0;

    return -M_models[flag]->getNeumannCondition(time, rate);
}

double
BCManager::
getOutflowNeumannJacobian(const double& time, const double& flag, const double& rate)
{
    auto it = M_models.find(flag);
    if (it == M_models.end())
        return 0.0;

    return -M_models[flag]->getNeumannJacobian(time, rate);
}

void
BCManager::
postProcess()
{
    for (auto windkessel : M_models)
        windkessel.second->shiftSolutions();
}

shp<VECTOREPETRA>
BCManager::
computeBoundaryIndicator(shp<FESPACE> fespace, const int& flag)
{
    shp<VECTOREPETRA> boundaryIndicator(new VECTOREPETRA(fespace->map()));
    boundaryIndicator->zero();

    LifeV::BCFunctionBase oneFunction(fOne);

    shp<LifeV::BCHandler> bcs;
    bcs.reset(new LifeV::BCHandler);

    bcs->addBC("BC", flag, LifeV::Essential, LifeV::Full, oneFunction, 3);

    bcs->bcUpdate(*fespace->mesh(), fespace->feBd(), fespace->dof());

    bcManageRhs(*boundaryIndicator, *fespace->mesh(), fespace->dof(),
                *bcs, fespace->feBd(), 1.0, 0.0);

    return boundaryIndicator;
}

double
BCManager::
fZero(const double& t, const double& x, const double& y,
      const double& z, const unsigned int& i)
{
    return 0.0;
}

double
BCManager::
fOne(const double& t, const double& x, const double& y,
      const double& z, const unsigned int& i)
{
    return 1.0;
}

double
BCManager::
constantFunction(const double& t, const double& x, const double& y,
                 const double& z, const unsigned int& i, const double& K)
{
    return K;
}

LifeV::MatrixSmall<3,3>
BCManager::
computeRotationMatrix(Vector3D vec, const unsigned int& index) const
{
    if (std::abs(vec.norm() - 1.0) <= 1e-10)
        vec.normalize();

    if (index > 2)
        printlog(YELLOW, "[computeRotationMatrix] WARNING: "
                         "the second argument must be either 0, 1 or 2!"
                         "Setting it to 2 as default value",
                 this->M_data.getVerbose());

    Matrix3D mat;
    Vector3D vec1;
    Vector3D vec2;

    // initialize vec1 and vec2 as orthogonal to vec
    vec1[0] = vec[1]; vec1[1] = -vec[0]; vec1[2] = 0.0;
    vec2[0] = vec[2]; vec2[1] = 0.0; vec2[2] = -vec[0];

    Vector3D ones(1.0);
    if (std::abs(vec1.norm()) <= 1e-10)
        vec1 = ones - vec - vec2;
    else if (std::abs(vec2.norm()) <= 1e-10)
        vec2 = ones - vec - vec1;

    // Gram-Schmidt (simplified as vec1 and vec2 are already orthogonal to vec)
    vec1.normalize();
    vec2 -= vec1 * (vec2.dot(vec1));
    vec2.normalize();

    // assign the prescribed normalized vector to the prescribed column
    mat[index][0] = vec[0];
    mat[index][1] = vec[1];
    mat[index][2] = vec[2];

    if (index != 0) {
        mat[0][0] = vec1[0];
        mat[0][1] = vec1[1];
        mat[0][2] = vec1[2];
    }
    else {
        mat[1][0] = vec1[0];
        mat[1][1] = vec1[1];
        mat[1][2] = vec1[2];
    }

    if (index != 2) {
        mat[2][0] = vec2[0];
        mat[2][1] = vec2[1];
        mat[2][2] = vec2[2];
    }
    else {
        mat[1][0] = vec2[0];
        mat[1][1] = vec2[1];
        mat[1][2] = vec2[2];
    }

    return mat;
}

std::map<unsigned int, LifeV::MatrixSmall<3,3>>
BCManager::
computeRotationMatrices() const
{
    std::map<unsigned int, Matrix3D> rotationMatrices;

    GeometricFace inlet = this->M_treeNode->M_block->getInlet();
    rotationMatrices[inlet.M_ringFlag] = this->computeRotationMatrix(inlet.M_normal, 2);

    std::vector<GeometricFace> outlets = this->M_treeNode->M_block->getOutlets();
    for (GeometricFace outlet : outlets)
       rotationMatrices[outlet.M_ringFlag] = this->computeRotationMatrix(outlet.M_normal, 2);

    return rotationMatrices;
}

shp<VECTOREPETRA>
BCManager::
computeRingsIndicator(shp<FESPACE> fespace) const
{
    shp<LifeV::BCHandler> bcs;
    bcs.reset(new LifeV::BCHandler);

    shp<VECTOREPETRA> ringsIndicator(new VECTOREPETRA(fespace->map()));
    ringsIndicator->zero();

    auto inletFlag = std::bind(constantFunction,
                               std::placeholders::_1,
                               std::placeholders::_2,
                               std::placeholders::_3,
                               std::placeholders::_4,
                               std::placeholders::_5,
                               M_inletRingFlag);
    LifeV::BCFunctionBase inletFunction(inletFlag);

    auto outletFlag = std::bind(constantFunction,
                                std::placeholders::_1,
                                std::placeholders::_2,
                                std::placeholders::_3,
                                std::placeholders::_4,
                                std::placeholders::_5,
                                M_outletRingFlag);
    LifeV::BCFunctionBase outletFunction(outletFlag);

    if (M_treeNode->isInletNode())
        bcs->addBC("InletRing", M_inletRingFlag, LifeV::EssentialEdges,
                   LifeV::Full, inletFunction, 3);

    if (M_treeNode->isOutletNode())
        bcs->addBC("OutletRing", M_outletRingFlag, LifeV::EssentialEdges,
                   LifeV::Full, outletFunction, 3);

    bcs->bcUpdate(*fespace->mesh(), fespace->feBd(), fespace->dof());

    bcManageRhs(*ringsIndicator, *fespace->mesh(), fespace->dof(),
            *bcs, fespace->feBd(), 1.0, 0.0);

    return ringsIndicator;
}

void
BCManager::
computeGlobalRotationMatrix(shp<FESPACE> fespace)
{
    // do not build the matrix if already built
    if (M_globalRotationMatrix)
        return;

    // throw Exception is the fespace is null
    if (!fespace)
        throw new Exception("Impossible to build a rotation matrix if the"
                            " given FESpace is null!");

    // define a proper matrix
    M_globalRotationMatrix.reset(new MATRIXEPETRA(fespace->map(), 3));
    M_globalRotationMatrix->insertOneDiagonal();  // define as the identity

    // compute the necessary rotation matrices
    std::map<unsigned int, Matrix3D> rotationMatrices = this->computeRotationMatrices();

    // defining a vector that identifies the ring DOFs
    shp<VECTOREPETRA> ringsIndicator = this->computeRingsIndicator(fespace);

    // initializing a temporary 3x3 matrix
    std::vector<double*> values(3);
    for (int n = 0; n < 3; n++)
    {
        values[n] = new double[3];
        for (int m = 0; m < 3; m++)
            values[n][m] = 0.0;
    }

    // initializing indices, rows and cols indices containers
    int indices[3];
    std::vector<int> rows;
    std::vector<int> cols;

    // Obtaining the ID of the elements
    const Epetra_Map epetraMap = *(fespace->mapPtr()->map(LifeV::Unique));
    unsigned int NumMyElements = epetraMap.NumMyElements();
    std::vector<int> MyGlobalElements(NumMyElements);
    epetraMap.MyGlobalElements(&MyGlobalElements[0]);

    unsigned int id;
    unsigned int flag = 999;

    unsigned int nDOF = fespace->dof().numTotalDof();

    for (unsigned int i = 0; i < NumMyElements/3; ++i)
    {
        id = MyGlobalElements[i];

        if (std::abs((*ringsIndicator)[id] - M_inletRingFlag) <= 1e-10)
            flag = M_inletRingFlag;
        else if (std::abs((*ringsIndicator)[id] - M_outletRingFlag) <= 1e-10)
            flag = M_outletRingFlag;
        else
            flag = 999;

        if (flag != 999)
        {
            indices[0] = id;
            indices[1] = id + nDOF;
            indices[2] = id + 2 * nDOF;

            cols.clear();
            cols.push_back (indices[0]);
            cols.push_back (indices[1]);
            cols.push_back (indices[2]);

            rows.clear();
            rows.push_back (indices[0]);
            rows.push_back (indices[1]);
            rows.push_back (indices[2]);

            // -1 because we added one to the diagonal
            M_globalRotationMatrix->addToCoefficient(indices[0], indices[0], -1.0);
            M_globalRotationMatrix->addToCoefficient(indices[1], indices[1], -1.0);
            M_globalRotationMatrix->addToCoefficient(indices[2], indices[2], -1.0);

            for (unsigned int r=0; r<3; r++) {
                for (unsigned int c=0; c<3; c++) {
                    values[r][c] = rotationMatrices[flag][r][c];
                }
            }

            M_globalRotationMatrix->addToCoefficients(3, 3, cols, rows, &values[0]);
        }
    }

    for (unsigned int n = 0; n < 3; n++)
        delete[] values[n];

    /*auto domainMap = M_globalRotationMatrix->domainMapPtr();
    auto rangeMap = M_globalRotationMatrix->rangeMapPtr();
    M_globalRotationMatrix->globalAssemble(domainMap, rangeMap);*/

    M_globalRotationMatrix->globalAssemble();
}

void
BCManager::
shiftToNormalTangentialCoordSystem(shp<MATRIXEPETRA> mat, shp<VECTOREPETRA> vec,
                                   shp<FESPACE> fespace)
{
    if ((fespace) && ((M_treeNode->isInletNode()) || (M_treeNode->isOutletNode())))
        this->computeGlobalRotationMatrix(fespace);

    unsigned int nDOF = fespace->dof().numTotalDof();

    if (mat)
    {
        unsigned int nRows = mat->matrixPtr()->NumGlobalRows();
        unsigned int nCols = mat->matrixPtr()->NumGlobalCols();

        // tmp1 = R*mat
        shp<MATRIXEPETRA> tmp1;
        if (nRows == 3 * nDOF)
        {
            tmp1.reset(new MATRIXEPETRA(mat->map(), mat->meanNumEntries()));
            M_globalRotationMatrix->multiply(false, *mat, false, *tmp1);
        }
        else
        {
            tmp1.reset(new MATRIXEPETRA(*mat));
        }

        // mat = tmp1*R^T"
        shp<MATRIXEPETRA> tmp2;
        if (nCols == 3 * nDOF)
        {
            tmp2.reset(new MATRIXEPETRA(mat->map(), mat->meanNumEntries()));
            tmp1->multiply(false, *(M_globalRotationMatrix->transpose()), false, *tmp2);

        }
        else
        {
            tmp2.reset(new MATRIXEPETRA(*tmp1));
        }

        mat->swapCrsMatrix(*tmp2);
    }

    if (vec)
    {
        unsigned int dim = vec->size();

        // vec = R*vec
        if (dim == 3 * nDOF)
        {
            shp<VECTOREPETRA> tmp;
            tmp.reset(new VECTOREPETRA(*vec));
            M_globalRotationMatrix->multiply(false, *tmp, *vec);
        }
    }
}

void
BCManager::
shiftToCartesianCoordSystem(shp<MATRIXEPETRA> mat, shp<VECTOREPETRA> vec,
                            shp<FESPACE> fespace)
{
    if ((fespace) && ((M_treeNode->isInletNode()) || (M_treeNode->isOutletNode())))
        this->computeGlobalRotationMatrix(fespace);

    unsigned int nDOF = fespace->dof().numTotalDof();

    if (mat)
    {
        unsigned int nRows = mat->matrixPtr()->NumGlobalRows();
        unsigned int nCols = mat->matrixPtr()->NumGlobalCols();

        // tmp1 = R^T*mat
        shp<MATRIXEPETRA> tmp1;
        if (nRows == 3 * nDOF)
        {
            tmp1.reset(new MATRIXEPETRA(mat->map(), mat->meanNumEntries()));
            M_globalRotationMatrix->transpose()->multiply(false, *mat, false, *tmp1);
        }
        else
        {
            tmp1.reset(new MATRIXEPETRA(*mat));
        }

        // mat = tmp1*R"
        shp<MATRIXEPETRA> tmp2;
        if (nCols == 3 * nDOF)
        {
            tmp2.reset(new MATRIXEPETRA(mat->map(), mat->meanNumEntries()));
            tmp1->multiply(false, *M_globalRotationMatrix, false, *tmp2);
        }
        else
        {
            tmp2.reset(new MATRIXEPETRA(*tmp1));
        }

        mat->swapCrsMatrix(*tmp2);
    }

    if (vec)
    {
        unsigned int dim = vec->size();

        // vec = R^T*vec
        if (dim == 3 * nDOF)
        {
            shp<VECTOREPETRA> tmp;
            tmp.reset(new VECTOREPETRA(*vec));
            M_globalRotationMatrix->transpose()->multiply(false, *tmp, *vec);
        }
    }
}


}
