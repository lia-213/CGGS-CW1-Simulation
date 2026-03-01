#ifndef SCENE_HEADER_FILE
#define SCENE_HEADER_FILE

#include <vector>
#include <fstream>
#include "ccd.h"
#include "volInt.h"
#include "auxfunctions.h"
#include "readMESH.h"
#include "mesh.h"
#include "constraints.h"

using namespace Eigen;
using namespace std;

// This class contains the entire scene operations, and the engine time loop.
class Scene {
public:
    double currTime;
    vector<Mesh> meshes;
    vector<Constraint> constraints;
    Mesh groundMesh;

    // Mostly for visualization
    MatrixXi allF, constEdges;
    MatrixXd currV, currConstVertices;

    // Adding objects. You do not need to update this generally.
    void add_mesh(const MatrixXd& V, const MatrixXi& F, const MatrixXi& T, const double density, const bool isFixed, const RowVector3d& COM, const RowVector4d& orientation) {
        Mesh m(V, F, T, density, isFixed, COM, orientation);
        meshes.push_back(m);

        MatrixXi newAllF(allF.rows() + F.rows(), 3);
        newAllF << allF, (F.array() + currV.rows()).matrix();
        allF = newAllF;
        MatrixXd newCurrV(currV.rows() + V.rows(), 3);
        newCurrV << currV, m.currV;
        currV = newCurrV;
    }

    /*********************************************************************
     This function handles a collision between objects m1 and m2 by assigning impulses.
     *********************************************************************/
    void handle_collision(Mesh& m1, Mesh& m2, const double& depth, const RowVector3d& contactNormal, const RowVector3d& penPosition, const double CRCoeff) {

        // n points from m1 into m2 (CCD convention). Ensure d > 0.
        double d = depth;
        RowVector3d n = contactNormal;
        if (d < 0) {
            d = -d;
            n = -n;
        }

        double invM1 = m1.isFixed ? 0.0 : m1.totalInvMass;
        double invM2 = m2.isFixed ? 0.0 : m2.totalInvMass;
        double sumInvMass = invM1 + invM2;

        if (sumInvMass == 0.0) return;

        // penPosition was adjusted in is_collide: pos -= depth*n/2
        // Recover the original CCD midpoint contact position
        RowVector3d contactPoint = penPosition + d * n / 2.0;

        // Compute moment arms from original COMs to the contact point
        RowVector3d r1 = contactPoint - m1.COM;
        RowVector3d r2 = contactPoint - m2.COM;

        // Velocity of each body at the contact point
        RowVector3d v1 = m1.comVelocity + m1.angVelocity.cross(r1);
        RowVector3d v2 = m2.comVelocity + m2.angVelocity.cross(r2);

        // n points m1->m2; vRelNormal > 0 means bodies approaching
        double vRelNormal = (v1 - v2).dot(n);

        // Resolve interpenetration: push m1 back along -n, m2 forward along +n
        double w1 = invM1 / sumInvMass;
        double w2 = invM2 / sumInvMass;
        m1.COM -= w1 * d * n;
        m2.COM += w2 * d * n;

        Matrix3d invIT1 = m1.isFixed ? Matrix3d::Zero() : m1.get_curr_inv_IT();
        Matrix3d invIT2 = m2.isFixed ? Matrix3d::Zero() : m2.get_curr_inv_IT();

        RowVector3d r1CrossN = r1.cross(n);
        RowVector3d r2CrossN = r2.cross(n);

        double rot1 = r1CrossN.dot((invIT1 * r1CrossN.transpose()).transpose());
        double rot2 = r2CrossN.dot((invIT2 * r2CrossN.transpose()).transpose());

        double denominator = sumInvMass + rot1 + rot2;

        // Standard collision impulse   
        double j = (1.0 + CRCoeff) * vRelNormal / denominator;
        RowVector3d impulse = j * n;

        // m1 pushed back along -n, m2 pushed forward along +n
        m1.comVelocity -= invM1 * impulse;
        m2.comVelocity += invM2 * impulse;

        m1.angVelocity -= (invIT1 * r1.cross(impulse).transpose()).transpose();
        m2.angVelocity += (invIT2 * r2.cross(impulse).transpose()).transpose();
    }

    /*********************************************************************
     This function handles a single time step.
     *********************************************************************/
    void update_scene(double timeStep, double CRCoeff, int maxIterations, double tolerance) {

        // 1. Integrate
        for (int i = 0; i < meshes.size(); i++)
            meshes[i].integrate(timeStep);

        // 2. Detect and Handle Mesh-Mesh Collisions
        double depth;
        RowVector3d contactNormal, penPosition;
        for (int i = 0; i < meshes.size(); i++) {
            for (int j = i + 1; j < meshes.size(); j++) {
                if (meshes[i].is_collide(meshes[j], depth, contactNormal, penPosition)) {
                    handle_collision(meshes[i], meshes[j], depth, contactNormal, penPosition, CRCoeff);
                }
            }
        }

        // 3. Detect and Handle Ground Collisions
        for (int i = 0; i < meshes.size(); i++) {
            int minyIndex;
            double minY = meshes[i].currV.col(1).minCoeff(&minyIndex);
            if (minY <= 0.0) {
                RowVector3d groundContactPoint = meshes[i].currV.row(minyIndex);
                groundContactPoint(1) = 0.0;
                handle_collision(meshes[i], groundMesh, minY, { 0.0, 1.0, 0.0 }, groundContactPoint, CRCoeff);
            }
        }

        // 4. Resolve Constraints
        int currIteration = 0;
        int zeroStreak = 0;
        int currConstIndex = 0;
        if (constraints.size() > 0) {
            while ((zeroStreak < constraints.size()) && (currIteration * constraints.size() < maxIterations)) {
                Constraint& currConstraint = constraints[currConstIndex];

                RowVector3d origConstPos1 = meshes[currConstraint.m1].origV.row(currConstraint.v1);
                RowVector3d origConstPos2 = meshes[currConstraint.m2].origV.row(currConstraint.v2);

                RowVector3d currConstPos1 = QRot(origConstPos1, meshes[currConstraint.m1].orientation) + meshes[currConstraint.m1].COM;
                RowVector3d currConstPos2 = QRot(origConstPos2, meshes[currConstraint.m2].orientation) + meshes[currConstraint.m2].COM;

                MatrixXd currCOMPositions(2, 3); currCOMPositions << meshes[currConstraint.m1].COM, meshes[currConstraint.m2].COM;
                MatrixXd currConstPositions(2, 3); currConstPositions << currConstPos1, currConstPos2;

                MatrixXd correctedCOMPositions;
                bool positionWasValid = currConstraint.resolve_position_constraint(currCOMPositions, currConstPositions, correctedCOMPositions, tolerance);

                if (positionWasValid) {
                    zeroStreak++;
                }
                else {
                    zeroStreak = 0;
                    meshes[currConstraint.m1].COM = correctedCOMPositions.row(0);
                    meshes[currConstraint.m2].COM = correctedCOMPositions.row(1);

                    currConstPos1 = QRot(origConstPos1, meshes[currConstraint.m1].orientation) + meshes[currConstraint.m1].COM;
                    currConstPos2 = QRot(origConstPos2, meshes[currConstraint.m2].orientation) + meshes[currConstraint.m2].COM;
                    currCOMPositions << meshes[currConstraint.m1].COM, meshes[currConstraint.m2].COM;
                    currConstPositions << currConstPos1, currConstPos2;

                    MatrixXd currCOMVelocities(2, 3); currCOMVelocities << meshes[currConstraint.m1].comVelocity, meshes[currConstraint.m2].comVelocity;
                    MatrixXd currAngVelocities(2, 3); currAngVelocities << meshes[currConstraint.m1].angVelocity, meshes[currConstraint.m2].angVelocity;

                    Matrix3d invInertiaTensor1 = meshes[currConstraint.m1].get_curr_inv_IT();
                    Matrix3d invInertiaTensor2 = meshes[currConstraint.m2].get_curr_inv_IT();
                    MatrixXd correctedCOMVelocities, correctedAngVelocities;

                    bool velocityWasValid = currConstraint.resolve_velocity_constraint(currCOMPositions, currConstPositions, currCOMVelocities, currAngVelocities, invInertiaTensor1, invInertiaTensor2, correctedCOMVelocities, correctedAngVelocities, tolerance);

                    if (!velocityWasValid) {
                        meshes[currConstraint.m1].comVelocity = correctedCOMVelocities.row(0);
                        meshes[currConstraint.m2].comVelocity = correctedCOMVelocities.row(1);
                        meshes[currConstraint.m1].angVelocity = correctedAngVelocities.row(0);
                        meshes[currConstraint.m2].angVelocity = correctedAngVelocities.row(1);
                    }
                }
                currIteration++;
                currConstIndex = (currConstIndex + 1) % (constraints.size());
            }
        }

        if (constraints.size() > 0 && currIteration * constraints.size() >= maxIterations)
            cout << "Constraint resolution reached maxIterations!" << endl;

        currTime += timeStep;

        // 5. Final Visual Update
        for (int i = 0; i < meshes.size(); i++) {
            for (int j = 0; j < meshes[i].currV.rows(); j++) {
                meshes[i].currV.row(j) << QRot(meshes[i].origV.row(j), meshes[i].orientation) + meshes[i].COM;
            }
        }

        int currVOffset = 0;
        for (int i = 0; i < meshes.size(); i++) {
            currV.block(currVOffset, 0, meshes[i].currV.rows(), 3) = meshes[i].currV;
            currVOffset += meshes[i].currV.rows();
        }
        for (int i = 0; i < constraints.size(); i += 2) {
            currConstVertices.row(i) = meshes[constraints[i].m1].currV.row(constraints[i].v1);
            currConstVertices.row(i + 1) = meshes[constraints[i].m2].currV.row(constraints[i].v2);
        }
    }

    // Loading scene (standard logic)
    bool load_scene(const std::string sceneFileName, const std::string constraintFileName) {
        ifstream sceneFileHandle;
        sceneFileHandle.open(DATA_PATH "/" + sceneFileName);
        if (!sceneFileHandle.is_open()) return false;

        int numofObjects;
        currTime = 0;
        sceneFileHandle >> numofObjects;
        for (int i = 0; i < numofObjects; i++) {
            MatrixXi objT, objF; MatrixXd objV;
            std::string MESHFileName; bool isFixed; double density;
            RowVector3d userCOM; RowVector4d userOrientation;
            sceneFileHandle >> MESHFileName >> density >> isFixed >> userCOM(0) >> userCOM(1) >> userCOM(2) >> userOrientation(0) >> userOrientation(1) >> userOrientation(2) >> userOrientation(3);
            userOrientation.normalize();
            readMESH(DATA_PATH "/" + MESHFileName, objV, objF, objT);
            MatrixXi tempF(objF.rows(), 3);
            tempF << objF.col(2), objF.col(1), objF.col(0);
            objF = tempF;
            add_mesh(objV, objF, objT, density, isFixed, userCOM, userOrientation);
        }
        groundMesh = Mesh(MatrixXd(0, 3), MatrixXi(0, 3), MatrixXi(0, 4), 0.0, true, RowVector3d::Zero(), RowVector4d::Zero());

        ifstream constraintFileHandle;
        constraintFileHandle.open(DATA_PATH "/" + constraintFileName);
        if (!constraintFileHandle.is_open()) return false;
        int numofConstraints;
        constraintFileHandle >> numofConstraints;
        currConstVertices.resize(numofConstraints * 2, 3);
        constEdges.resize(numofConstraints, 2);
        for (int i = 0; i < numofConstraints; i++) {
            int attachM1, attachM2, attachV1, attachV2;
            double lowerBound, upperBound;
            constraintFileHandle >> attachM1 >> attachV1 >> attachM2 >> attachV2 >> lowerBound >> upperBound;
            double initDist = (meshes[attachM1].currV.row(attachV1) - meshes[attachM2].currV.row(attachV2)).norm();
            double invMass1 = meshes[attachM1].isFixed ? 0.0 : meshes[attachM1].totalInvMass;
            double invMass2 = meshes[attachM2].isFixed ? 0.0 : meshes[attachM2].totalInvMass;
            constraints.push_back(Constraint(DISTANCE, INEQUALITY, false, attachM1, attachV1, attachM2, attachV2, invMass1, invMass2, RowVector3d::Zero(), lowerBound * initDist, 0.0));
            constraints.push_back(Constraint(DISTANCE, INEQUALITY, true, attachM1, attachV1, attachM2, attachV2, invMass1, invMass2, RowVector3d::Zero(), upperBound * initDist, 0.0));
            currConstVertices.row(2 * i) = meshes[attachM1].currV.row(attachV1);
            currConstVertices.row(2 * i + 1) = meshes[attachM2].currV.row(attachV2);
            constEdges.row(i) << 2 * i, 2 * i + 1;
        }
        return true;
    }

    Scene() { allF.resize(0, 3); currV.resize(0, 3); }
    ~Scene() {}
};

#endif