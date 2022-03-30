#include "Bypass.hpp"

namespace RedMA
{

    Bypass::
    Bypass(EPETRACOMM comm, std::string name, bool verbose, bool boundary_layer, bool randomizable) :
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
        computeStenosisCenter();
        computeStenosisOuterNormal();

        M_identity3D(0,0) = 1;
        M_identity3D(0,1) = 0;
        M_identity3D(0,2) = 0;
        M_identity3D(1,0) = 0;
        M_identity3D(1,1) = 1;
        M_identity3D(1,2) = 0;
        M_identity3D(2,0) = 0;
        M_identity3D(2,1) = 0;
        M_identity3D(2,2) = 1;

        setDistorsionMatrix();
    }

    void
    Bypass::
    resetInletOutlets()
    {
        GeometricFace inlet1(M_inletCenterRef1, M_inletNormalRef1, M_inletRadiusRef1, 2, 20);
        GeometricFace inlet2(M_inletCenterRef2, M_inletNormalRef2, M_inletRadiusRef2, 3, 30);
        GeometricFace outlet(M_outletCenterRef, M_outletNormalRef, M_outletRadiusRef, 4, 40);

        M_inlets.clear();
        M_inlets.push_back(inlet1);
        M_inlets.push_back(inlet2);
        M_outlets.clear();
        M_outlets.push_back(outlet);
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
    computeStenosisCenter()
    {
        M_stenosisCenter[0] = -7.8933;
        M_stenosisCenter[1] = 2.5276;
        M_stenosisCenter[2] = 42.4082;
    }

    void
    Bypass::
    computeStenosisOuterNormal()
    {
        M_stenosisOuterNormal[0] = 0.32901267;
        M_stenosisOuterNormal[1] = 0.94276643;
        M_stenosisOuterNormal[2] = -0.05425531;
    }

    void
    Bypass::
    setDistorsionMatrix()
    {
        // setting the eigenvector matrix

        // x eigenvector for the geometry
        M_Eigenvector1[0] = 0.85784321;
        M_Eigenvector1[1] = -0.32731597;
        M_Eigenvector1[2] = -0.30261594;

        // y eigenvector for the geometry
        M_Eigenvector2 = M_stenosisOuterNormal;

        // z eigenvector for the geometry
        M_Eigenvector3[0] = 0.30307373;
        M_Eigenvector3[1] = -0.05905023;
        M_Eigenvector3[2] = 0.95113583;

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

        diagMatrix(0,0) = 2;
        diagMatrix(0,1) = 0;
        diagMatrix(0,2) = 0;
        diagMatrix(1,0) = 0;
        diagMatrix(1,1) = 2;
        diagMatrix(1,2) = 0;
        diagMatrix(2,0) = 0;
        diagMatrix(2,1) = 0;
        diagMatrix(2,2) = 1;

        // distorsion matrix for the stenosis norm
        M_distorsionMatrix = eigenMatrix * diagMatrix * eigenMatrix.inverse();
    }

    double
    Bypass::
    stenosisDeformation(const double &z) {
        return (std::abs(z) < 1) * (0.5 + 0.5 * std::sin(std::atan(1) * 4 * (z + 0.5)));
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
                shp <Transformer> transformer, bool transformMesh) {

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

                shp <LifeV::BCHandler> bcs(new LifeV::BCHandler);
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
                nAffineDeformer.deformMesh(*transformer);
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
                nAffineDeformer.deformMesh(*transformer);
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

        addStenosis(M_parametersHandler["stenosis_amplitude"],
                    M_parametersHandler["stenosis_width"],
                    transformer, transformMesh);
        
        M_mesh->check(1, true);

        bend(M_parametersHandler["in1_alphax"],
             M_parametersHandler["in1_alphay"],
             M_parametersHandler["in1_alphaz"],
             M_parametersHandler["in2_alphax"],
             M_parametersHandler["in2_alphay"],
             M_parametersHandler["in2_alphaz"], transformer, transformMesh);

        printlog(MAGENTA, "done\n", M_verbose);
    }

}