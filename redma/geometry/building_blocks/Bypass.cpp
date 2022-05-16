#include "Bypass.hpp"

namespace RedMA
{

    Bypass::
    Bypass(EPETRACOMM comm, std::string name, bool verbose, bool boundary_layer, bool isBifurcation, unsigned int activeStenosis,
           bool randomizable):
            BuildingBlock(comm, "coarse", verbose)
    {
        M_name = name;
        M_datafileName = "data_mesh";
        if (!(std::strcmp(M_refinement.c_str(), "coarse")))
        {
            if (boundary_layer)
                M_meshName = "others/bypass_BL.mesh";
            else
                M_meshName = "others/bypass_coarse_fluid.mesh";
        }
        else
            throw new Exception("No refined meshed for the bypass geometry are available! "
                                "Please set the refinement to coarse.");

        // center of inlet1 (reference configuration)
        M_inletCenterRef1[0] = -7.06006787;
        M_inletCenterRef1[1] = 1.63241487;
        M_inletCenterRef1[2] = 45.36727996;

        // center of inlet2 (reference configuration)
        M_inletCenterRef2[0] = -8.51731429;
        M_inletCenterRef2[1] = 2.71285607;
        M_inletCenterRef2[2] = 45.23909379;

        // center of outlet (reference configuration)
        M_outletCenterRef[0] = -8.49107192;
        M_outletCenterRef[1] = 2.23961129;
        M_outletCenterRef[2] = 40.33716354;

        // normal of inlet1 (reference configuration)
        M_inletNormalRef1[0] = 0.12364816;
        M_inletNormalRef1[1] = -0.21329954;
        M_inletNormalRef1[2] = 0.96913076;

        // outlet of inlet2 (reference configuration)
        M_inletNormalRef2[0] = -0.13280704;
        M_inletNormalRef2[1] = 0.30002537;
        M_inletNormalRef2[2] = 0.94464124;

        // outlet of outlet (reference configuration)
        M_outletNormalRef[0] = 0.27216076;
        M_outletNormalRef[1] = 0.1600235;
        M_outletNormalRef[2] = 0.94885246;

        // inlets and outlet radia (reference configuration)
        M_inletRadiusRef1 = 0.33609344929489127;
        M_inletRadiusRef2 = 0.2661432021360946;
        M_outletRadiusRef = 0.24480087734232522;

        M_isBifurcation = isBifurcation;

        M_wallFlag = 200;

        resetInletOutlets();

        const double maxAngle = 0.4;
        const double maxAmplitude = 0.4;
        const double maxWidth = 0.4;

        M_parametersHandler.registerParameter("in1_alphax", 0.0, -maxAngle,
                                              maxAngle, randomizable);
        M_parametersHandler.registerParameter("in1_alphay", 0.0, -maxAngle,
                                              maxAngle, randomizable);
        M_parametersHandler.registerParameter("in1_alphaz", 0.0, -maxAngle,
                                              maxAngle, randomizable);
        M_parametersHandler.registerParameter("in2_alphax", 0.0, -maxAngle,
                                              maxAngle, randomizable);
        M_parametersHandler.registerParameter("in2_alphay", 0.0, -maxAngle,
                                              maxAngle, randomizable);
        M_parametersHandler.registerParameter("in2_alphaz", 0.0, -maxAngle,
                                              maxAngle, randomizable);
        M_parametersHandler.registerParameter("stenosis_amplitude", 0.0, 0,
                                              maxAmplitude, randomizable);
        M_parametersHandler.registerParameter("stenosis_width", 0.0, 0,
                                              maxWidth, randomizable);

        computeCenter();


        M_identity3D(0,0) = 1;
        M_identity3D(0,1) = 0;
        M_identity3D(0,2) = 0;
        M_identity3D(1,0) = 0;
        M_identity3D(1,1) = 1;
        M_identity3D(1,2) = 0;
        M_identity3D(2,0) = 0;
        M_identity3D(2,1) = 0;
        M_identity3D(2,2) = 1;

        setStenosisAttributes();
        setActiveStenosis(activeStenosis);
        setDistorsionMatrix();
    }

    void
    Bypass::
    resetInletOutlets()
    {

        M_inlets.clear();
        M_outlets.clear();

        if (M_isBifurcation)
        {
            GeometricFace inlet1(M_inletCenterRef1, -1.0 * M_inletNormalRef1, M_inletRadiusRef1, 2, 20);
            GeometricFace inlet2(M_inletCenterRef2, -1.0 * M_inletNormalRef2, M_inletRadiusRef2, 3, 30);
            GeometricFace outlet(M_outletCenterRef, -1.0 * M_outletNormalRef, M_outletRadiusRef, 4, 40);

            M_inlets.push_back(outlet);
            M_outlets.push_back(inlet1);
            M_outlets.push_back(inlet2);
        }
        else
        {
            GeometricFace inlet1(M_inletCenterRef1, M_inletNormalRef1, M_inletRadiusRef1, 2, 20);
            GeometricFace inlet2(M_inletCenterRef2, M_inletNormalRef2, M_inletRadiusRef2, 3, 30);
            GeometricFace outlet(M_outletCenterRef, M_outletNormalRef, M_outletRadiusRef, 4, 40);

            M_inlets.push_back(inlet1);
            M_inlets.push_back(inlet2);
            M_outlets.push_back(outlet);
        }
    }

    void
    Bypass::
    computeCenter()
    {
        /*double angleInlet1 = std::atan(M_outlet1NormalRef[1] /
                M_outlet1NormalRef[2]);
        double zCenter = M_outlet1NormalRef[2] -
                M_outlet1CenterRef[1] * std::tan(angleInlet1);
        M_center[0] = 0;
        M_center[1] = 0;
        M_center[2] = zCenter;*/

        M_center[0] = -7.95993;
        M_center[1] = 2.18022;
        M_center[2] = 42.7164;
    }

    void
    Bypass::
    setStenosisAttributes()
    {

        // inlet
        Vector3D inletStenosisCenter;
        inletStenosisCenter[0] = -8.29497;
        inletStenosisCenter[1] = 2.11155;
        inletStenosisCenter[2] = 41.1408;
        Vector3D inletStenosisNormal;
        inletStenosisNormal[0] = -0.0883044;
        inletStenosisNormal[1] = -0.995009;
        inletStenosisNormal[2] = 0.0464609;
        Vector3D inletStenosisFirstEigenvector;
        inletStenosisFirstEigenvector[0] = -0.9427618;
        inletStenosisFirstEigenvector[1] = 0.0934536;
        inletStenosisFirstEigenvector[2] = 0.02974148;
        Vector3D inletStenosisThirdEigenvector;
        inletStenosisThirdEigenvector[0] = 0.0456323;
        inletStenosisThirdEigenvector[1] = 0.04070732;
        inletStenosisThirdEigenvector[2] = 0.98426454;
        std::map<std::string, Vector3D> mapInletStenosis;
        mapInletStenosis.insert(std::pair<std::string, Vector3D> ("center", inletStenosisCenter));
        mapInletStenosis.insert(std::pair<std::string, Vector3D> ("normal", inletStenosisNormal));
        mapInletStenosis.insert(std::pair<std::string, Vector3D> ("eigenX", inletStenosisFirstEigenvector));
        mapInletStenosis.insert(std::pair<std::string, Vector3D> ("eigenY", inletStenosisNormal));
        mapInletStenosis.insert(std::pair<std::string, Vector3D> ("eigenZ", inletStenosisThirdEigenvector));
        M_stenosisAttributes.insert(std::pair<unsigned int, std::map<std::string, Vector3D>> (0, mapInletStenosis));


        // old stenosis before branches
        Vector3D StenosisCenter;
        StenosisCenter[0] = -7.8933;
        StenosisCenter[1] = 2.5276;
        StenosisCenter[2] = 42.4082;
        Vector3D StenosisNormal;
        StenosisNormal[0] = 0.32901267;
        StenosisNormal[1] = 0.94276643;
        StenosisNormal[2] = -0.05425531;
        Vector3D StenosisFirstEigenvector;
        StenosisFirstEigenvector[0] = 0.85784321;
        StenosisFirstEigenvector[1] = -0.32731597;
        StenosisFirstEigenvector[2] = -0.30261594;
        Vector3D StenosisThirdEigenvector;
        StenosisThirdEigenvector[0] = 0.30307373;
        StenosisThirdEigenvector[1] = -0.05905023;
        StenosisThirdEigenvector[2] = 0.95113583;
        std::map<std::string, Vector3D> mapStenosis;
        mapStenosis.insert(std::pair<std::string, Vector3D> ("center", StenosisCenter));
        mapStenosis.insert(std::pair<std::string, Vector3D> ("normal", StenosisNormal));
        mapStenosis.insert(std::pair<std::string, Vector3D> ("eigenX", StenosisFirstEigenvector));
        mapStenosis.insert(std::pair<std::string, Vector3D> ("eigenY", StenosisNormal));
        mapStenosis.insert(std::pair<std::string, Vector3D> ("eigenZ", StenosisThirdEigenvector));
        M_stenosisAttributes.insert(std::pair<unsigned int, std::map<std::string, Vector3D>> (1, mapStenosis));


        //stenosis at first outlet
        Vector3D newStenosisCenter;
        newStenosisCenter[0] = -7.30198;
        newStenosisCenter[1] = 2.22133;
        newStenosisCenter[2] = 44.35341;
        Vector3D newStenosisNormal;
        newStenosisNormal[0] = -0.132643;
        newStenosisNormal[1] = 0.947123;
        newStenosisNormal[2] = 0.292172;
        Vector3D newStenosisFirstEigenvector;
        newStenosisFirstEigenvector[0] = 0.9681090;
        newStenosisFirstEigenvector[1] = 0.1948699;
        newStenosisFirstEigenvector[2] = -0.2111276;
        Vector3D newStenosisThirdEigenvector;
        newStenosisThirdEigenvector[0] = 0.25590652;
        newStenosisThirdEigenvector[1] = -0.2638001;
        newStenosisThirdEigenvector[2] = 0.94567632;
        std::map<std::string, Vector3D> mapNewStenosis;
        mapNewStenosis.insert(std::pair<std::string, Vector3D> ("center", newStenosisCenter));
        mapNewStenosis.insert(std::pair<std::string, Vector3D> ("normal", newStenosisNormal));
        mapNewStenosis.insert(std::pair<std::string, Vector3D> ("eigenX", newStenosisFirstEigenvector));
        mapNewStenosis.insert(std::pair<std::string, Vector3D> ("eigenY", newStenosisNormal));
        mapNewStenosis.insert(std::pair<std::string, Vector3D> ("eigenZ", newStenosisThirdEigenvector));
        M_stenosisAttributes.insert(std::pair<unsigned int, std::map<std::string, Vector3D>> (2, mapNewStenosis));


        Vector3D secondNewStenosisCenter;
        secondNewStenosisCenter[0] = -8.99104;
        secondNewStenosisCenter[1] = 2.45104;
        secondNewStenosisCenter[2] = 44.71721;
        Vector3D secondNewStenosisNormal;
        secondNewStenosisNormal[0] = 0.992063;
        secondNewStenosisNormal[1] = 0.125465;
        secondNewStenosisNormal[2] = -0.008358;
        Vector3D secondNewStenosisFirstEigenvector;
        secondNewStenosisFirstEigenvector[0] = 0.115863;
        secondNewStenosisFirstEigenvector[1] = -0.884866;
        secondNewStenosisFirstEigenvector[2] = 0.385514;
        Vector3D secondNewStenosisThirdEigenvector;
        secondNewStenosisThirdEigenvector[0] = -0.0526753;
        secondNewStenosisThirdEigenvector[1] = 0.395785;
        secondNewStenosisThirdEigenvector[2] = 0.909554;
        std::map<std::string, Vector3D> mapSecondNewStenosis;
        mapSecondNewStenosis.insert(std::pair<std::string, Vector3D> ("center", secondNewStenosisCenter));
        mapSecondNewStenosis.insert(std::pair<std::string, Vector3D> ("normal", secondNewStenosisNormal));
        mapSecondNewStenosis.insert(std::pair<std::string, Vector3D> ("eigenX", secondNewStenosisNormal));
        mapSecondNewStenosis.insert(std::pair<std::string, Vector3D> ("eigenY", secondNewStenosisFirstEigenvector));
        mapSecondNewStenosis.insert(std::pair<std::string, Vector3D> ("eigenZ", secondNewStenosisThirdEigenvector));
        M_stenosisAttributes.insert(std::pair<unsigned int, std::map<std::string, Vector3D>> (3, mapSecondNewStenosis));
    }


    void
    Bypass::
    setActiveStenosis(unsigned int i)
    {
       if (i > 3)
           throw new Exception("There are 4 stenosis available but a number bigger than 3 was given as input...");
       std::string msg = "Selecting stenosis " + std::to_string(i) + "\n";
       printlog(GREEN, msg);
       M_stenosisCenter = M_stenosisAttributes[i]["center"];
       M_stenosisOuterNormal = M_stenosisAttributes[i]["normal"];
       M_Eigenvector1 = M_stenosisAttributes[i]["eigenX"];
       M_Eigenvector2 = M_stenosisAttributes[i]["eigenY"];
       M_Eigenvector3 = M_stenosisAttributes[i]["eigenZ"];
       if (i == 0)
           M_diameterAtStenosis = 0.4912600;
       else if (i == 1)
           M_diameterAtStenosis = 0.6325955;
       else if (i == 2)
           M_diameterAtStenosis = 0.6585998;
       else if (i == 3)
           M_diameterAtStenosis = 0.5109350;
    }

    void
    Bypass::
    setDistorsionMatrix()
    {
        // setting the eigenvector matrix
        Matrix3D eigenMatrix;

        eigenMatrix(0,0) = M_Eigenvector1[0];
        eigenMatrix(0,1) = M_Eigenvector2[0];
        eigenMatrix(0,2) = M_Eigenvector3[0];
        eigenMatrix(1,0) = M_Eigenvector1[1];
        eigenMatrix(1,1) = M_Eigenvector2[1];
        eigenMatrix(1,2) = M_Eigenvector3[1];
        eigenMatrix(2,0) = M_Eigenvector1[2];
        eigenMatrix(2,1) = M_Eigenvector2[2];
        eigenMatrix(2,2) = M_Eigenvector3[2];

        // setting the eigenvalues matrix
        Matrix3D diagMatrix;

        diagMatrix(0,0) = 1;
        diagMatrix(0,1) = 0;
        diagMatrix(0,2) = 0;
        diagMatrix(1,0) = 0;
        diagMatrix(1,1) = 1;
        diagMatrix(1,2) = 0;
        diagMatrix(2,0) = 0;
        diagMatrix(2,1) = 0;
        diagMatrix(2,2) = 0.5;

        // distorsion matrix for the stenosis norm
        M_distorsionMatrix = eigenMatrix * diagMatrix * eigenMatrix.inverse();
    }

    double
    Bypass::
    stenosisDeformation(const double &z) {
        return (std::abs(z) < 1) * (std::exp(1 / (pow(z, 2)-1)));
    }

    double
    Bypass::
    stenosisBC(const double &t, const double &x,
               const double &y, const double &z,
               const LifeV::ID &i, const double& amplitude, const double& width,
               const Vector3D& stenosisCentre, const Vector3D& stenosisNormal,
               const Matrix3D& distorsionMatrix) {

        Vector3D point(x, y, z);
        Vector3D distorsionVector = distorsionMatrix * (point-stenosisCentre);
        double deformationNorm = (point-stenosisCentre).dot(distorsionVector);
        return - 1.0 * amplitude * stenosisDeformation(deformationNorm / width) * stenosisNormal[i];
    }

    void
    Bypass::
    addStenosis(const double &amplitude, const double &width,
                shp <Transformer> transformer, bool transformMesh)
    {
        using namespace std::placeholders;

        if (amplitude > 0 && width > 0) {

            std::string msg = std::string("[") + M_name + " BuildingBlock]";
            msg = msg + " introducing a stenosis with amplitude = (" + std::to_string(amplitude)
                  + ") and width = (" + std::to_string(width) + " ) \n";
            printlog(GREEN, msg, M_verbose);

            auto fooWall = std::bind(stenosisBC,
                                     std::placeholders::_1,
                                     std::placeholders::_2,
                                     std::placeholders::_3,
                                     std::placeholders::_4,
                                     std::placeholders::_5,
                                     amplitude, width, M_stenosisCenter, M_stenosisOuterNormal, M_distorsionMatrix);

            if (transformMesh) {
                NonAffineDeformer nAffineDeformer(M_mesh, M_comm, M_verbose);

                LifeV::BCFunctionBase zeroFunction(BuildingBlock::fZero);
                LifeV::BCFunctionBase stenosisFunction(fooWall);

                shp<LifeV::BCHandler> bcs(new LifeV::BCHandler);
                bcs->addBC("Inlet1", 2, LifeV::Essential, LifeV::Full,
                           zeroFunction, 3);
                bcs->addBC("Inlet2", 3, LifeV::Essential, LifeV::Full,
                           zeroFunction, 3);
                bcs->addBC("Outlet", 4, LifeV::Essential, LifeV::Full,
                           zeroFunction, 3);
                bcs->addBC("Wall", 200, LifeV::Essential, LifeV::Full,
                           stenosisFunction, 3);

                CoutRedirecter ct;
                ct.redirect();
                nAffineDeformer.applyBCs(bcs);
                std::string xmlFilename = M_datafile("geometric_structure/xmldeformer",
                                                     "SolverParamList.xml");
                nAffineDeformer.setXMLsolver(xmlFilename);
                M_displacement = nAffineDeformer.solveSystem("Ifpack");
                nAffineDeformer.deformMeshComposite(*transformer, M_displacement);
                printlog(CYAN, ct.restore(), M_verbose);
            }

        }

    }

    void
    Bypass::
    bend(const double& in1_alphax, const double& in1_alphay, const double& in1_alphaz,
         const double& in2_alphax, const double& in2_alphay, const double& in2_alphaz,
         shp<Transformer> transformer, bool transformMesh)
    {
        if ((std::abs(in1_alphax) > 0 || std::abs(in1_alphay) > 0 || std::abs(in1_alphaz) > 0) ||
            (std::abs(in2_alphax) > 0 || std::abs(in2_alphay) > 0 || std::abs(in2_alphaz) > 0))
        {
            std::string msg = std::string("[") + M_name + " BuildingBlock]";
            msg = msg + " bending with angles = (" + std::to_string(in1_alphax)
                  + ", " + std::to_string(in1_alphay) +
                  ", " + std::to_string(in1_alphaz) + ") at outlet 1, and (" +
                  std::to_string(in2_alphax) + ", " +
                  std::to_string(in2_alphay) + ", " +
                  std::to_string(in2_alphaz) + ")"
                  + " at outlet 2" + "\n";
            printlog(GREEN, msg, M_verbose);

            using namespace std::placeholders;

            Vector3D rotationCenter = M_center;

            // handle rotation for inlet 1
            Matrix3D rotationMatrix =
                    computeRotationMatrix(0, in1_alphax) *
                    computeRotationMatrix(1, in1_alphay) *
                    computeRotationMatrix(2, in1_alphaz);

            Vector3D rotatedCenterInlet1;
            Vector3D rotatedNormalInlet1;
            rotateGeometricFace(M_inlets[0], rotatedCenterInlet1, rotatedNormalInlet1,
                                rotationMatrix, rotationCenter);

            auto fooInlet1 = std::bind(inletMapFunction,
                                       std::placeholders::_1,
                                       std::placeholders::_2,
                                       std::placeholders::_3,
                                       std::placeholders::_4,
                                       std::placeholders::_5,
                                       M_inlets[0], rotatedCenterInlet1,
                                       rotationMatrix);

            M_inlets[0].M_center = rotatedCenterInlet1;
            M_inlets[0].M_normal = rotatedNormalInlet1;

            // handle rotation for inlet 2
            rotationMatrix = computeRotationMatrix(0, in2_alphax) *
                             computeRotationMatrix(1, in2_alphay) *
                             computeRotationMatrix(2, in2_alphaz);

            Vector3D rotatedCenterInlet2;
            Vector3D rotatedNormalInlet2;
            rotateGeometricFace(M_inlets[1], rotatedCenterInlet2, rotatedNormalInlet2,
                                rotationMatrix, rotationCenter);

            auto fooInlet2 = std::bind(inletMapFunction,
                                       std::placeholders::_1,
                                       std::placeholders::_2,
                                       std::placeholders::_3,
                                       std::placeholders::_4,
                                       std::placeholders::_5,
                                       M_inlets[1], rotatedCenterInlet2,
                                       rotationMatrix);

            M_inlets[1].M_center = rotatedCenterInlet2;
            M_inlets[1].M_normal = rotatedNormalInlet2;

            if (transformMesh)
            {
                NonAffineDeformer nAffineDeformer(M_mesh, M_comm, M_verbose);

                LifeV::BCFunctionBase zeroFunction(BuildingBlock::fZero);
                LifeV::BCFunctionBase inletFunction1(fooInlet1);
                LifeV::BCFunctionBase inletFunction2(fooInlet2);

                shp<LifeV::BCHandler> bcs(new LifeV::BCHandler);
                bcs->addBC("Inlet1", 2, LifeV::Essential, LifeV::Full,
                           inletFunction1, 3);
                bcs->addBC("Inlet2", 3, LifeV::Essential, LifeV::Full,
                           inletFunction2, 3);
                bcs->addBC("Outlet", 4, LifeV::Essential, LifeV::Full,
                           zeroFunction, 3);

                CoutRedirecter ct;
                ct.redirect();
                nAffineDeformer.applyBCs(bcs);
                std::string xmlFilename = M_datafile("geometric_structure/xmldeformer",
                                                     "SolverParamList.xml");
                nAffineDeformer.setXMLsolver(xmlFilename);
                printlog(CYAN, ct.restore(), M_verbose);
            }
        }
    }

    double
    Bypass::
    inletMapFunction(const double& t, const double& x,
                     const double& y, const double& z,
                     const LifeV::ID& i, const GeometricFace& targetFace,
                     const Vector3D& desiredCenter, const Matrix3D& rotationMatrix)
    {
        Vector3D cur(x, y, z);
        Vector3D diff = cur - targetFace.M_center;

        Vector3D modifiedDif = rotationMatrix * diff;
        Vector3D newPoint = desiredCenter + modifiedDif;

        return newPoint[i] - cur[i];
    }

    void
    Bypass::
    rotateGeometricFace(const GeometricFace& face, Vector3D& rotatedCenter,
                        Vector3D& rotatedNormal, const Matrix3D& rotationMatrix,
                        const Vector3D& rotationCenter)
    {
        rotatedNormal = rotationMatrix * face.M_normal;
        rotatedCenter = rotationMatrix * (face.M_center - rotationCenter) +
                        rotationCenter;
    }


    std::string
    Bypass::
    getOptionalParameter(unsigned int index)
    {
    }

    void
    Bypass::
    applyNonAffineTransformation(bool transformMesh)
    {
        printlog(MAGENTA, "[" + M_name + " BuildingBlock] applying non affine transformation ...\n", M_verbose);

        shp<LifeV::MeshUtility::MeshTransformer<MESH> > transformer;
        if (transformMesh)
            transformer.reset(new Transformer(*M_mesh));

        addStenosis(M_parametersHandler["stenosis_amplitude"] * M_diameterAtStenosis,
                    M_parametersHandler["stenosis_width"],
                    transformer, transformMesh);

        transformer->savePoints();
        printlog(MAGENTA, "done\n", M_verbose);
    }

}