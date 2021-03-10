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

    M_ringConstraint = data("bc_conditions/ring_constraint", "normal");

    if (M_treeNode->isOutletNode())
        this->parseOutflowNeumannData();
    if (M_treeNode->isInletNode())
        this->checkInflowLaw();

    M_inletFlag = treeNode->M_block->getInlet().M_flag;
    for(auto out_face : treeNode->M_block->getOutlets())
        // these are all the outlets, even the ones with children in the block structure!
        M_outletFlags.push_back(out_face.M_flag);
    for(auto out_face : treeNode->getOutlets())
        // these are the "true" outlets, i.e. faces without children!
        M_trueOutletFlags.push_back(out_face.M_flag);
    M_wallFlag = treeNode->M_block->getWallFlag();
    M_inletRingFlag = treeNode->M_block->getInlet().M_ringFlag;
    for (auto out_face : treeNode->M_block->getOutlets())
        // these are all the outlets, even the ones with children in the block structure!
        M_outletRingFlags.push_back(out_face.M_ringFlag);
    for (auto out_face : treeNode->getOutlets())
        // these are the "true" outlets, i.e. faces without children!
        M_trueOutletRingFlags.push_back(out_face.M_ringFlag);
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

    auto foo = std::bind(this->poiseuilleInflow,
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
applyOutflowNeumannBCs(shp<LifeV::BCHandler> bcs, const bool& zeroFlag) const
{
    std::function<double(double)> outflowLaw;
    if (zeroFlag)
        outflowLaw = std::function<double(double)>([](double t){return 0.0;});
    else
        // TODO: set an outflow function is non-homo Neumann outflow BCs wish to be imposed!
        outflowLaw = std::function<double(double)>([](double t){return 0.0;});

    auto outflowBoundaryCondition = std::bind(neumannInflow,
                                              std::placeholders::_1,
                                              std::placeholders::_2,
                                              std::placeholders::_3,
                                              std::placeholders::_4,
                                              std::placeholders::_5,
                                              outflowLaw);

    LifeV::BCFunctionBase outflowFunction(outflowBoundaryCondition);

    for (unsigned int flag : M_outletFlags)
        bcs->addBC("OutletHomoNeumann", flag, LifeV::Natural, LifeV::Normal,
                   outflowFunction);
}

void
BCManager::
addInletBC(shp<LifeV::BCHandler> bcs, const bool& ringOnly, const bool& zeroFlag) const
{
    if (M_treeNode->isInletNode())
    {
        if (!std::strcmp(M_inletBCType.c_str(), "dirichlet"))
            // imposing an inflow velocity profile
            this->applyInflowDirichletBCs(bcs, zeroFlag);

        else if (!std::strcmp(M_inletBCType.c_str(), "neumann"))
            // imposing an inflow pressure --> actually this is done at RHS so this is not called!
            this->applyInflowNeumannBCs(bcs, zeroFlag);

        else
            throw new Exception("Selected inflow BCs not implemented!");

         /*LifeV::BCFunctionBase zeroFunction(fZero);

        if (!ringOnly) {
            bcs->addBC("Wall", M_wallFlag, LifeV::Essential,
                       LifeV::Full, zeroFunction, 3);
        }
        else {
            bcs->addBC("InletRing", M_inletRingFlag, LifeV::EssentialEdges,
                       LifeV::Full, zeroFunction, 3);

            std::vector <LifeV::ID> compz (1);
            compz[0] = 2;
            bcs->addBC("InletRingZ", M_inletRingFlag, LifeV::EssentialEdges,
                       LifeV::Component, zeroFunction, compz);
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
        shp<LifeV::BCHandler> bcsRing = createBCHandler0DirichletRing();
        bcsRing->bcUpdate(*fespace->mesh(), fespace->feBd(), fespace->dof());

        bool shiftReferenceSystem = (!index);

        if (shiftReferenceSystem && (!std::strcmp(M_ringConstraint.c_str(), "normal")))
            // shift to (t1,t2,n) reference system
            this->shiftToNormalTangentialCoordSystem(nullptr, curVec, fespace);

        // apply proper ring BCs (either normal or full)
        bcManageRhs(*curVec, *fespace->mesh(), fespace->dof(),
                    *bcsRing, fespace->feBd(), 0.0, 0.0);

        if (shiftReferenceSystem && (!std::strcmp(M_ringConstraint.c_str(), "normal")))
            // shift back to (x,y,z) reference system
            this->shiftToCartesianCoordSystem(nullptr, curVec, fespace);
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

    unsigned int nRows = input.nRows();
    unsigned int nCols = input.nCols();
    for (unsigned int j = 0; j < nCols; j++)
    {
        if (!input.block(index, j)->isZero())
        {
            std::cout << "Applying bcs to " << index << " " << j << std::endl << std::flush;
            shp<MATRIXEPETRA> curMatrix = spcast<MATRIXEPETRA>(input.block(index, j)->data());
            auto domainMap = curMatrix->domainMapPtr();
            auto rangeMap = curMatrix->rangeMapPtr();

            bcManageMatrix(*curMatrix, *fespace->mesh(),
                           fespace->dof(), *bcs, fespace->feBd(),
                           (j == index) * diagCoefficient, 0.0);

            if (ringOnly)
            {
                shp<LifeV::BCHandler> bcsRing = createBCHandler0DirichletRing();
                bcsRing->bcUpdate(*fespace->mesh(), fespace->feBd(), fespace->dof());

                // the rotation is needed only on 'velocity x velocity' matrices!
                bool shiftReferenceSystem = (nRows==2) && (nCols==2) && (!index) && (!j);

                if (shiftReferenceSystem && (!std::strcmp(M_ringConstraint.c_str(), "normal")))
                    // shift to (t1,t2,n) reference system
                    this->shiftToNormalTangentialCoordSystem(curMatrix, nullptr, fespace);

                // apply BCs
                bcManageMatrix(*curMatrix, *fespace->mesh(),
                               fespace->dof(), *bcsRing, fespace->feBd(),
                               (j == index) * diagCoefficient, 0.0);

                if (shiftReferenceSystem && (!std::strcmp(M_ringConstraint.c_str(), "normal")))
                    // shift back to in (x,y,z) reference system
                    this->shiftToCartesianCoordSystem(curMatrix, nullptr, fespace);
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
        shp<LifeV::BCHandler> bcsRing = createBCHandler0DirichletRing();
        bcsRing->bcUpdate(*fespace->mesh(), fespace->feBd(), fespace->dof());

        bool shiftReferenceSystem = (!index);

        if (shiftReferenceSystem && (!std::strcmp(M_ringConstraint.c_str(), "normal")))
            // shift to (t1,t2,n) reference system
            this->shiftToNormalTangentialCoordSystem(nullptr, curVec, fespace);

        // apply proper ring BCs (either normal or full)
        bcManageRhs(*curVec, *fespace->mesh(), fespace->dof(),
                    *bcsRing, fespace->feBd(), 0.0, 0.0);

        if (shiftReferenceSystem && (!std::strcmp(M_ringConstraint.c_str(), "normal")))
            // shift back to (x,y,z) reference system
            this->shiftToCartesianCoordSystem(nullptr, curVec, fespace);
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

    if (!std::strcmp(M_ringConstraint.c_str(), "normal")) {
        if (M_treeNode->isInletNode())
        {
            bcs->addBC("InletRingZ", M_inletRingFlag, LifeV::EssentialEdges,
                       LifeV::Component, zeroFunction, compz);
        }

        if (M_treeNode->isOutletNode())
        {
            for(const unsigned int& ringFlag : M_trueOutletRingFlags)
                bcs->addBC("OutletRingZ", ringFlag, LifeV::EssentialEdges,
                           LifeV::Component, zeroFunction, compz);
        }
    }
    else if (!std::strcmp(M_ringConstraint.c_str(), "full"))
    {
        if (M_treeNode->isInletNode())
        {
            bcs->addBC("InletRing", M_inletRingFlag, LifeV::EssentialEdges,
                       LifeV::Full, zeroFunction, 3);
        }

        if (M_treeNode->isOutletNode())
        {
            for(const unsigned int& ringFlag : M_trueOutletRingFlags)
                bcs->addBC("OutletRing", ringFlag, LifeV::EssentialEdges,
                           LifeV::Full, zeroFunction, 3);
        }
    }
    else
    {
        throw new Exception("Unrecognized type of ring constraint!"
                            "Recognized types: {normal, full}");
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

std::vector<unsigned int>
BCManager::
getWallFlags(const bool& withRings) const
{
    std::vector<unsigned int> flags;

    flags.push_back(M_wallFlag);
    if (withRings)
    {
        flags.push_back(M_inletRingFlag);
        flags.insert(std::end(flags), std::begin(M_outletRingFlags), std::end(M_outletRingFlags));
    }

    return flags;
}

shp<VECTOREPETRA>
BCManager::
computeBoundaryIndicator(shp<FESPACE> fespace, const std::vector<unsigned int> flags) const
{
    shp<VECTOREPETRA> boundaryIndicator(new VECTOREPETRA(fespace->map()));
    boundaryIndicator->zero();

    LifeV::BCFunctionBase oneFunction(fOne);

    shp<LifeV::BCHandler> bcs;
    bcs.reset(new LifeV::BCHandler);

    for (unsigned int flag : flags)
    {
        if ((flag == M_inletFlag) || (flag == M_wallFlag) ||
            (std::find(M_outletFlags.begin(), M_outletFlags.end(), flag) != M_outletFlags.end()))
            bcs->addBC("BC", flag, LifeV::Essential, LifeV::Full, oneFunction, 3);
        else if ((flag == M_inletRingFlag) ||
                 (std::find(M_outletRingFlags.begin(), M_outletRingFlags.end(), flag) != M_outletRingFlags.end()))
            bcs->addBC("BC", flag, LifeV::EssentialEdges, LifeV::Full, oneFunction, 3);
    }

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

    std::pair<Vector3D, Vector3D> tangents = this->computeTangentVersors(vec);

    // assign the prescribed normalized vector to the prescribed column
    mat[0][index] = vec[0];
    mat[1][index] = vec[1];
    mat[2][index] = vec[2];

    if (index != 0) {
        mat[0][0] = tangents.first[0];
        mat[1][0] = tangents.first[1];
        mat[2][0] = tangents.first[2];
    }
    else {
        mat[0][1] = tangents.first[0];
        mat[1][1] = tangents.first[1];
        mat[2][1] = tangents.first[2];
    }

    if (index != 2) {
        mat[0][2] = tangents.second[0];
        mat[1][2] = tangents.second[1];
        mat[2][2] = tangents.second[2];
    }
    else {
        mat[0][1] = tangents.second[0];
        mat[1][1] = tangents.second[1];
        mat[2][1] = tangents.second[2];
    }

    return mat;
}

std::pair<LifeV::VectorSmall<3>, LifeV::VectorSmall<3>>
BCManager::
computeTangentVersors(const Vector3D& normal) const
{
    Vector3D t1;
    Vector3D t2;
    std::pair<Vector3D, Vector3D> res;

    double nx = normal[0];
    double ny = normal[1];
    double nz = normal[2];

    double nx2 = std::sqrt(normal[1] * normal[1] + normal[2] * normal[2]);
    double ny2 = std::sqrt(normal[0] * normal[0] + normal[2] * normal[2]);
    double nz2 = std::sqrt(normal[0] * normal[0] + normal[1] * normal[1]);

    if ((nx2 >= ny2) && (nx2 >= nz2))
    {
        //We create t1
        t1[0] = 0;
        t1[1] = nz / nx2;
        t1[2] = -ny / nx2;

        //We create t2
        t2[0] = -nx2;
        t2[1] = nx * ny / nx2;
        t2[2] = nx * nz / nx2;
    }
    else if ((ny2 >= nx2) && (ny2 >= nz2))
    {
        //We create t1
        t1[0] = -nz / ny2;
        t1[1] = 0;
        t1[2] = nx / ny2;

        //We create t2
        t2[0] = nx * ny / ny2;
        t2[1] = -ny2;
        t2[2] = ny * nz / ny2;
    }
    else
    {
        //We create t1
        t1[0] = ny / nz2;
        t1[1] = -nx / nz2;
        t1[2] = 0;

        //We create t2
        t2[0] = nx * nz / nz2;
        t2[1] = ny * nz / nz2;
        t2[2] = -nz2;
    }

    res.first = t1;
    res.second = t2;

    return res;
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
    typedef std::function<double(const double&, const double&, const double&, const double&, const unsigned int&)> FUN;

    shp<LifeV::BCHandler> bcs;
    bcs.reset(new LifeV::BCHandler);

    shp<VECTOREPETRA> ringsIndicator(new VECTOREPETRA(fespace->map()));
    ringsIndicator->zero();

    if (M_treeNode->isInletNode())
    {
        FUN inletFlag = std::bind(constantFunction,
                                  std::placeholders::_1,
                                  std::placeholders::_2,
                                  std::placeholders::_3,
                                  std::placeholders::_4,
                                  std::placeholders::_5,
                                  M_inletRingFlag);
        LifeV::BCFunctionBase inletFunction(inletFlag);

        bcs->addBC("InletRing", M_inletRingFlag, LifeV::EssentialEdges,
                   LifeV::Full, inletFunction, 3);
    }

    std::string baseBCName = "OutletRing_";
    std::string BCName;

    if (M_treeNode->isOutletNode())
    {
        FUN outletFlag;
        LifeV::BCFunctionBase outletFunction;

        for(unsigned int ringFlag : M_trueOutletRingFlags)
        {
            outletFlag = std::bind(constantFunction,
                                   std::placeholders::_1,
                                   std::placeholders::_2,
                                   std::placeholders::_3,
                                   std::placeholders::_4,
                                   std::placeholders::_5,
                                   ringFlag);
            outletFunction.setFunction(outletFlag);

            BCName = baseBCName + std::to_string(ringFlag);
            bcs->addBC(BCName, ringFlag, LifeV::EssentialEdges,
                       LifeV::Full, outletFunction, 3);
        }
    }

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
        else
        {
            auto it = std::find(M_trueOutletRingFlags.begin(), M_trueOutletRingFlags.end(), (*ringsIndicator)[id]);
            if (it != M_trueOutletRingFlags.end())
                flag = *it;
            else
                flag = 999;
        }

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
